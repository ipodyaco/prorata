#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
import time
import itertools, copy, math
import subprocess
import shlex


from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import MCS
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys



from pubchempy import *

def parse_options(argv):

    opts, args = getopt.getopt(argv[1:], "hi:o:n:s:",
                                    ["help",
                             	     "input-file",
	                			     "output-file",
                                     "decoy-number",
                                     "similarity-threshold"])

    output_filename = ""
    input_filename  = ""
    iDecoyNumber    = 0
    dSimilarity_threshold     = -1

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-file, -o output-file -n decoy-number -s similarity-threshold"
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-file"):
            output_filename = value
        if option in ("-n", "decoy-number"):
            iDecoyNumber = int(value)
        if option in ("-s", "similarity-threshold"):
            dSimilarity_threshold = float(value)

    if (input_filename == "") :
        print "Please specify -i"
        sys.exit(1)
    if (iDecoyNumber  <= 0) :
        print "Please specify -n"
        sys.exit(1)
    if (dSimilarity_threshold < 0):
        print "Please sepcify -s"
        sys.exit(1)

    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_pubchemdecoy" + inputFileNameExt
    

    return (input_filename, output_filename, iDecoyNumber, dSimilarity_threshold)

def compareInchi(sInchi_first, sInchi_second):
    
    Mol_first  = Chem.MolFromInchi(sInchi_first)
    try:
        Mol_second = Chem.MolFromInchi(sInchi_second)
        Mols_list  = [Mol_first, Mol_second]
        fps = [MACCSkeys.GenMACCSKeys(x) for x in Mols_list]
        dSimilarity = DataStructs.FingerprintSimilarity(fps[0],fps[1])
    except:
        print "bad inchi", sInchi_second
        dSimilarity = 2
    return dSimilarity

def GenerateDecoy(inchi_ori, iDecoyNumber, dSimilarity_threshold) :
    inchi_ele_list = inchi_ori.split("/")
    formula_ori    = inchi_ele_list[1]
    #print inchi_ori
    all_same_formula_cids_list = get_cids(formula_ori, 'formula')
    #all_same_formula_compounds_list = get_compounds(formula_ori, 'formula')
    decoy_list = []
    for each_cid in all_same_formula_cids_list :
        current_compound_list = get_compounds(each_cid, 'cid')
        if (len(current_compound_list) != 1) :
            print current_compound_list
            sys.exit(1)
        else :
            current_decoy_inchi = current_compound_list[0].inchi.encode('ascii','ignore')
        if (compareInchi(inchi_ori, current_decoy_inchi) < dSimilarity_threshold) :
            decoy_list.append(current_decoy_inchi)
            if (len(decoy_list) == iDecoyNumber) :
                break
    if (len(decoy_list) < iDecoyNumber) :
        print "Only find", len(decoy_list), "decoys for", formula_ori
    return decoy_list

def GenerateNewDB(input_filename, output_filename, iDecoyNumber, dSimilarity_threshold) :
    input_file = open(input_filename)
    output_file= open(output_filename, "w")
    #temp_filename = os.path.dirname(output_filename)+os.sep+"omgtemp.sdf"

    for each_line in input_file:
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        output_file.write(each_line+"\n")
        if (each_line.startswith( "#")) :
            continue
        compound_info_list  = each_line.split("\t")
        sInternalID_ori = compound_info_list[0] 
        inchi_ori  = compound_info_list[1]
        #decoy_compound_list = GenerateDecoy(inchi_ori, iDecoyNumber, temp_filename)
        decoy_compound_list = GenerateDecoy(inchi_ori, iDecoyNumber, dSimilarity_threshold)
        decoy_count = 0
        for each_decoy_inchi in decoy_compound_list :

            sInternalID_decoy = "Decoy_"+str(decoy_count)+"_"+sInternalID_ori
            decoy_count += 1
            output_file.write(sInternalID_decoy+"\t"+each_decoy_inchi)
            for i in range(2, len(compound_info_list)):
                output_file.write("\t"+compound_info_list[i])
            output_file.write("\n")
    input_file.close()
    output_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        input_filename ,  output_filename, iDecoyNumber, dSimilarity_threshold = parse_options(argv)  

    GenerateNewDB(input_filename, output_filename, iDecoyNumber, dSimilarity_threshold)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()

