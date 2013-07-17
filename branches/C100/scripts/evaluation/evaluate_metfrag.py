#!/usr/bin/python

import sys, getopt, warnings, os, re
import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hw:o:r:",
                                    ["help",
                                     "working-dir",
				                     "output-dir",
                                     "realhit-filename"])


    # Default working dir and config file
    working_dir      = ""
    output_dir       = ""
    realhit_filename = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w working-dir -o output-dir -r realhit_filename"
            sys.exit(0)
        if option in ("-w", "--working-dir"):
            working_dir  = value
            if (working_dir[-1] != os.sep) :
                working_dir += os.sep
        if option in ("-o", "--output-dir"):
            output_dir   = value
            if (output_dir[-1] != os.sep) :
                output_dir += os.sep
        if option in ("-r", "--realhit-filename"):
            realhit_filename  = value

    if ( working_dir == ""  ) or (realhit_filename == "")  :
        print "please specify working dir and realhit filename"
        sys.exit(1)

    if (output_dir == "") :
        output_dir = working_dir
    
    sdf_filename_list = get_file_list_with_ext(working_dir, ".sdf")

    return [working_dir, output_dir, realhit_filename]

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

def extract_realhits(realhit_filename) :
    realhit_list = []
    realhit_file = open(realhit_filename)
    for each_line in realhit_file :
        each_line = each_line.strip()
        if (each_line == ""):
            continue
        if (each_line.startswith("*\t")) :
            scan_info = each_line.split("\t")
            scan_id   = scan_info[1]
            scan_inchi= scan_info[4]
            realhit_list.append([scan_id, scan_inchi])
    realhit_file.close()
    return realhit_list

def evalue_one_prediction(scan_id, scan_inchi, working_dir, output_dir) :
    sdf_result_filename = working_dir+"results_"+scan_id+".sdf"
    organized_result_filename = output_dir+"results_"+scan_id+".txt"
    sdf_mol_list = Chem.SDMolSupplier(sdf_result_filename)
    all_results_list = []
    organized_result_file = open(organized_result_filename, "w")
    for current_sdf_mol in sdf_mol_list :
        dScore = float(current_sdf_mol.GetProp("Score"))
        sInchi = current_sdf_mol.GetProp("Inchi")
        sUniqueID = current_sdf_mol.GetProp("UniqueID")
        #print  sUniqueID, sInchi,  dScore
        all_results_list.append([dScore, sUniqueID, sInchi])
    all_results_list.sort(key=lambda e:e[0], reverse=True)
    realhit_rank = 0
    output_details = ""
    for i in range(len(all_results_list)) :
        dScore    = all_results_list[i][0]
        sUniqueID = all_results_list[i][1]
        sInchi    = all_results_list[i][2]
        if (i==0) :
            current_rank = 1
        elif (dScore < dPrevious_Score) :
            current_rank = i+1
        output_details += str(current_rank)+"\t"+sUniqueID+"\t"+str(dScore)+"\t"+sInchi+"\n"
        dPrevious_Score = dScore
        if (sInchi == scan_inchi) :
            realhit_rank = current_rank
            realhit_score= dScore
    if (realhit_rank == 0) :
        print "can't find the realhit", scan_id, scan_inchi
        realhit_rank = len(all_results_list) +1
        realhit_score= -100.0
    organized_result_file.write(">"+str(realhit_rank)+"\t"+scan_id+"\t"+str(realhit_score)+"\t"+scan_inchi+"\n")
    organized_result_file.write(output_details)
    organized_result_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [working_dir, output_dir, realhit_filename] = parse_options(argv)      
        realhit_list =  extract_realhits(realhit_filename) 
        for scan_id, scan_inchi in realhit_list :
            evalue_one_prediction(scan_id, scan_inchi, working_dir, output_dir)

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



