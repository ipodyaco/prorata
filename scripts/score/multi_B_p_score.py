#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math
from multiprocessing import Pool

import scoring_B

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hvbri:p:o:e:n:t:",
                                    ["help",
                                     "verbose",
                                     "break-ring",
                                     "ranksum",
                                     "input-filename",
                                     "compound-filename",
				                     "output-dir",
                                     "energy-filename",
                                     "process-number",
                                     "precursor-type"])

    # Default working dir and config file
    input_filename  = ""
    compound_filename = ""
    output_foldername   = ""
    energy_filename   = ""
    process_number    = 0
    bSpectrumDetails  = False
    bBreakRing        = False
    precursor_type    = 1
    bRankSum          = False
    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print " -i input_filename -p compound_filename -o folder_dir -e energy_filename -n process_number -v (verbose) -b (breakring) -t precursor-type -r (ranksum)"
            sys.exit(0)
        if option in ("-v", "--verbose") :
            bSpectrumDetails = True
        if option in ("-b", "--break-ring") :
            bBreakRing  = True
        if option in ("-i", "--input-filename"):
            input_filename  = value
        if option in ("-r", "--ranksum"):
            bRankSum = True
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
        if option in ("-t", "--precursor-type") :
            precursor_type = int(value)
            if ((precursor_type != 1) and (precursor_type != -1)) :
                print "-t precursor type must be 1 or -1"
                sys.exit(1)

    if ((input_filename == "") or (compound_filename == "") or (output_foldername == "") or (energy_filename == "") or (process_number == 0)) :
        print "please specify working dir, input file name, compound file, output dir, and energy file"
        sys.exit(1)
    
    return [input_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing, precursor_type, bRankSum]


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

def HandleOneRealHit (script_dir, compound_filename, output_filename, energy_filename, bSpectrumDetails, bBreakRing, precursor_type, bRankSum, precursor_mz, allPeaks_list, s_chemical_structure) :
#    print output_filename, s_chemical_structure    
    scoring_B.score_main(compound_filename, output_filename, energy_filename, bSpectrumDetails, bBreakRing, precursor_type, bRankSum, precursor_mz, allPeaks_list, s_chemical_structure)


def HandleAllRealHits(script_dir, all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing, precursor_type, bRankSum):
    dProtonMass = 1.007825 # proton mass
    mypool = Pool(processes=process_number)
    all_realhit_file = open(all_realhit_filename)

    for each_line in all_realhit_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith("*")) :
            allPeaks_list = []
            each_realhit_corename = each_line.split("\t")[1].strip()
            current_output_filename = output_foldername + each_realhit_corename  + ".output.txt"
            precursor_mz = float(each_line.split("\t")[3].strip())
            s_chemical_structure = each_line.split("\t")[4].strip()
            #print script_dir+os.sep, current_output_filename  
        if (each_line.startswith("+")) :
            current_mofz      = float(each_line[1:].split(" ")[0].strip())
            current_intensity = float(each_line[1:].split(" ")[-1].strip())
            allPeaks_list.append([current_mofz, current_intensity])
        if (each_line.startswith("//")) :            
            result = mypool.apply_async(HandleOneRealHit,(script_dir, compound_filename,current_output_filename,energy_filename, bSpectrumDetails, bBreakRing, precursor_type, bRankSum, precursor_mz, allPeaks_list, s_chemical_structure ))
            precursor_mz = ""
            s_chemical_structure = 0
    mypool.close()
    mypool.join()
#    if result.successful():
#        print 'successful'

    all_realhit_file.close()

def main(argv=None):

    #precursor_accuracy = 0.02
    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing, precursor_type, bRankSum] = parse_options(argv)
        script_path = sys.argv[0]
        script_dir  = os.path.dirname(os.path.abspath(script_path))
        HandleAllRealHits(script_dir, all_realhit_filename, compound_filename, output_foldername, energy_filename, process_number, bSpectrumDetails, bBreakRing, precursor_type, bRankSum)



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()


