#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time

from rdkit import Chem
from rdkit.Chem import Descriptors

def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:c:o:",
                                    ["help",
                                     "working-dir",
				                     "output-file",])


    # Default working dir and config file
    working_dir = "./"
    outputFileName = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory -o outputfile"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-file"):
            outputFileName = value

    SMILES_filename_list = get_file_list_with_ext(working_dir, ".txt")
    
    if (outputFileName == "") :
        outputFileName = working_dir + "output"

    return [SMILES_filename_list ,  outputFileName]

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
        die("Program exit!")

    return file_list

def parseMassbankFile(currentInchiFileName) :
    current_Inchi_file = open(currentInchiFileName)

    for each_line in current_Inchi_file :    
        each_line = each_line.strip()
        if each_line.startswith("CH$IUPAC:") :
            s_chemical_structure = each_line.split(" ")[1]
            print s_chemical_structure
            break
    current_Inchi_file.close()
    return s_chemical_structure

def BreakOneBond(current_mol, iBondsNum) :
    for i in range(iBondsNum) :
        current_editable_mol = Chem.EditableMol(current_mol)
        current_bond  = current_mol.GetBondWithIdx(i)
        idx_beginAtom = current_bond.GetBeginAtomIdx()
        idx_endAtom   = current_bond.GetEndAtomIdx()
        current_editable_mol.RemoveBond(idx_beginAtom, idx_endAtom)
        current_modified_mol = current_editable_mol.GetMol()
        current_fragments    = Chem.GetMolFrags(current_modified_mol, asMols=True)
        #current_fragments    = Chem.GetMolFrags(current_modified_mol)
        if (len(current_fragments) > 1) :
            for each_fragment in current_fragments :
                current_sInchi = Chem.MolToInchi(each_fragment)
                current_dMass  = Descriptors.ExactMolWt(each_fragment)
                current_sBondType = current_bond.GetBondType()
                print current_dMass, current_sBondType, current_sInchi
        #        print each_fragment
def ExhaustBonds(sInchiInfo) :
    current_mol = Chem.MolFromInchi(sInchiInfo)
    iBondsNum = current_mol.GetNumBonds()
    BreakOneBond(current_mol, iBondsNum)
    

def HandleInchi(currentInchiFileName, output_file) :
    output_file.write(">\t"+currentInchiFileName+"\n")
    sInchiInfo = parseMassbankFile(currentInchiFileName)
    ExhaustBonds(sInchiInfo)
    


def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [Inchi_filename_list ,  outputFileName] = parse_options(argv)  
    output_file = open(outputFileName, "w")
    for eachInchiFileName in Inchi_filename_list : 
        HandleInchi(eachInchiFileName, output_file)

    output_file.close()



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
