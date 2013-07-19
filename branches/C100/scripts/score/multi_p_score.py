#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math
from multiprocessing import Pool

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hvbw:i:p:o:e:n:",
                                    ["help",
                                     "verbose",
                                     "break-ring",
                                     "working-dir",
                                     "input-filename",
                                     "compound-filename",
				                     "output-dir",
                                     "energy-filename",
                                     "process-number"])

    # Default working dir and config file
    working_dir     = ""
    input_filename  = ""
    compound_filename = ""
    output_foldername   = ""
    energy_filename   = ""
    process_number    = 0
    bSpectrumDetails  = False
    bBreakRing        = False

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w working-dir -i input_filename -p compound_filename -o folder_dir -e energy_filename -n process_number -v verbose -b breakring"
            sys.exit(0)
        if option in ("-w", "--working-dir") :
            working_dir = value
            if (working_dir[-1] != os.sep) :
                working_dir += os.sep
        if option in ("-v", "--verbose") :
            bSpectrumDetails = True
        if option in ("-b", "--break-ring") :
            bBreakRing  = True
        if option in ("-i", "--input-filename"):
            input_filename  = value
        if option in ("-p", "--compound-filename"):
            compound_filename = value
        if option in ("-o", "--output-dir"):
            output_foldername = value
            if (output_foldername[-1:]!="/") and (output_foldername != ""):
                output_foldername += "/"
        if option in ("-e", "--energy-filename"):
            energy_filename   = value
        if option in ("-n", "--process-number"):
            process_number = int(value)

    if ((working_dir == "") or (input_filename == "") or (compound_filename == "") or (output_foldername == "") or (energy_filename == "") or (process_number == 0)) :
        print "please specify working dir, input file name, compound file, output dir, and energy file"
        sys.exit(1)
    
    return [working_dir, input_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing]


## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        sys.exit(1)

    return file_list

def HandleOneRealHit (script_dir, realhit_filename, compound_filename, output_filename, energy_filename, bSpectrumDetails, bBreakRing) :
    # change scoring_A.py to scoring.py will go to metfrag scoring function
    if (bSpectrumDetails) :
        vstr = " -v "
    else :
        vstr = ""
    if (bBreakRing) :
        bstr = " -b "
    else :
        bstr = ""
    command_str = "python "+script_dir+os.sep+"scoring_A.py -m "+realhit_filename+" -p "+compound_filename+" -o "+output_filename+" -e "+energy_filename + vstr + bstr
#    print command_str
    os.system(command_str)

def HandleAllRealHits(script_dir, working_dir, all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing):
   
    mypool = Pool(processes=process_number)
    all_realhit_file = open(all_realhit_filename)

    for each_line in all_realhit_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith("*")) :
            each_realhit_corename = each_line.split("\t")[1].strip()
            each_realhit_filename = working_dir + each_realhit_corename + ".txt"
            if not (os.path.isfile(each_realhit_filename)) :
                print "can't find", each_realhit_filename
                sys.exit(1)
            current_output_filename = output_foldername + os.path.basename(each_realhit_filename) + ".output.txt"
            #print script_dir+os.sep, current_output_filename     
            result = mypool.apply_async(HandleOneRealHit,(script_dir,each_realhit_filename,compound_filename,current_output_filename,energy_filename, bSpectrumDetails, bBreakRing))
    mypool.close()
    mypool.join()
#    if result.successful():
#        print 'successful'

    all_realhit_file.close()

def main(argv=None):

    precursor_accuracy = 0.02
    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [working_dir, all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing] = parse_options(argv)
        script_path = sys.argv[0]
        script_dir  = os.path.dirname(os.path.abspath(script_path))
        HandleAllRealHits(script_dir, working_dir, all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing)



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()


