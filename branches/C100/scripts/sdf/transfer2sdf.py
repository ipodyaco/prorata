#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                                     "input-filename",
				                     "output-filename"])

    # Default working dir and config file
    input_filename  = ""
    output_filename = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input_filename -o output_filename"
            sys.exit(0)
        if option in ("-i", "--input-filename"):
            input_filename  = value
        if option in ("-o", "--output-filename"):
            output_filename   = value

    if input_filename == "" :
        print "please specify input file name"
        sys.exit(1)
    if output_filename == "" :
        (input_root, input_ext) = os.path.splitext(input_filename)
        output_filename = input_root + ".sdf"
    
    return [input_filename, output_filename]

def transferFile(input_filename, output_filename) :

    input_file = open(input_filename)
    writer = Chem.SDWriter(output_filename)

    for each_line in input_file :
        each_line = each_line.strip()
        each_line_info_list = each_line.split("\t")
        sUniqueID   = each_line_info_list[0]
        sInchi_data = each_line_info_list[1]
        sDBLink     = each_line_info_list[3]
        current_mol = Chem.MolFromInchi(sInchi_data)
        current_mol.SetProp("UniqueID", sUniqueID)
        current_mol.SetProp("DBLink", sDBLink)
        writer.write(current_mol)

    writer.close()
    input_file.close()

def main(argv=None):
    if argv is None:
        argv = sys.argv
       		 # parse options
        [input_filename, output_filename] = parse_options(argv)  	
        transferFile(input_filename, output_filename)

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



