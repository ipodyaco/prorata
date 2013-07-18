#!/usr/bin/python

import sys, getopt, warnings, os, re
import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hr:o:",
                                    ["help",
                                     "realhit-filename",
				                     "output-filename"])


    # Default working dir and config file
    realhit_filename  = ""
    output_filename   = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-r realhit-filename -o output_filename"
            sys.exit(0)
        if option in ("-o", "--output-filename"):
            output_filename   = value
        if option in ("-r", "--realhit-filename"):
            realhit_filename = value

    if (output_filename == "") or (realhit_filename == "") :
        print "please specify realhit filename and output filename"
        sys.exit(1)
    
    return [realhit_filename,  output_filename]

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


def calculateMass(sInchi) :
    current_mol  = Chem.MolFromInchi(sInchi)
    dCalculatedMass =  Descriptors.ExactMolWt(current_mol)
    return dCalculatedMass

def extractMassBankFile(input_filename, output_file, all_inchi) :

    
    input_file = open(input_filename)

    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith("*")) :
            compound_info = each_line.split("\t")
            sAccession = compound_info[1]
            sPubChemId = compound_info[2]
            sInchi     = compound_info[4]
            sCompoundName = compound_info[5]
            dCalculatedMass = calculateMass(sInchi)
        if (sInchi not in all_inchi) :
            all_inchi.append(sInchi)
            output_file.write("MassBank_"+sAccession+"\t"+sInchi+"\t"+str(dCalculatedMass)+"\t"+sCompoundName)
            if (sPubChemId != "NA"):
                output_file.write("\t(PUBCHEM \""+sPubChemId+"\")\n")
            else :
                output_file.write("\t"+sPubChemId+"\n")

    input_file.close()

#    output_file.write("*\t"+sAccession+"\t"+sPubChemId+"\t"+str(dPrecursorEM)+"\t"+sInchi+"\n")
#    for each_peak in sPeak_list :
#        output_file.write("+"+each_peak[0]+" "+each_peak[1]+"\n")
#    output_file.write("//\n")

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [realhit_filename, output_filename] = parse_options(argv)  
        #badcase_dir = working_dir+"nopubchemid"
        #os.system("mkdir " + badcase_dir)
        output_file = open(output_filename, "w")
        all_inchi = [] 
        extractMassBankFile(realhit_filename, output_file, all_inchi)

        output_file.close()

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



