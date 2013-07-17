#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
import time
import itertools, copy, math
import subprocess
import shlex
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def parse_options(argv):

    opts, args = getopt.getopt(argv[1:], "hi:o:n:",
                                    ["help",
                             	     "input-file",
	                			     "output-dir",
                                     "decoy-number"])

    output_dir = ""
    input_filename  = ""
    iDecoyNumber    = 0

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-file, -o output-dir -n decoy-number"
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-dir"):
            output_dir = value
        if option in ("-n", "decoy-number"):
            iDecoyNumber = int(value)

    if ((input_filename == "")  or (output_dir == "")):
        print "Please specify -i and -o"
        sys.exit(1)
    if (iDecoyNumber  <= 0) :
        print "Please specify -n"
        sys.exit(1)

    return (input_filename, output_dir, iDecoyNumber)

def CallOMG(ori_formula, temp_filename) :
#    command_str = os.environ["JAVA_HOME"]+"/bin/java -jar PMG_1.0.jar "+ ori_formula +" -filter -p 1 -o "+temp_filename
   # sJAVA_Head  = "/home/yf1/software/java/1.6/jdk1.6.0_37"
    sJAVA_Head  = os.environ["JAVA_HOME"]
    command_str = "java -jar "+os.environ["OMG_HOME"]+"/OMG.jar -ec "+ ori_formula+" -o " + temp_filename
    myarg = shlex.split(command_str)
    #print myarg
    #os.system(command_str)
    proc = subprocess.Popen(myarg)
    while (True) :
        time.sleep(60)
        temp_file_size = os.path.getsize(temp_filename)
        if (temp_file_size > 1000000):
            pstatus = proc.poll()
            if pstatus is None:
                proc.kill()
            break
    proc.wait()

def GenerateDecoy_mp(sInternalID_ori, inchi_ori, iDecoyNumber, ori_line, output_dir, compound_info_list) :
    temp_filename = output_dir+os.sep+sInternalID_ori+".sdf"
    decoy_list = []
    inchi_list = [inchi_ori]
    ori_mol = Chem.MolFromInchi(inchi_ori)
    ori_formula = AllChem.CalcMolFormula(ori_mol)
    ori_fragments_list = Chem.GetMolFrags(ori_mol, asMols=True, sanitizeFrags=False)
    if (len(ori_fragments_list) == 1) :
        bCheckFragmentNumber = True
    else :
        bCheckFragmentNumber = False
    CallOMG(ori_formula, temp_filename)
    decoy_mol_list = Chem.SDMolSupplier(temp_filename)
    used_decoy_count = 0
    for each_decoy_mol in decoy_mol_list :
        if each_decoy_mol is None :
            continue
        current_inchi = Chem.MolToInchi(each_decoy_mol)
        if (current_inchi in inchi_list) :
            continue
        if (bCheckFragmentNumber) :
            decoy_fragments_list = Chem.GetMolFrags(each_decoy_mol, asMols=True, sanitizeFrags=False)
            if (len(decoy_fragments_list) != 1) :
                continue
        decoy_list.append(current_inchi)
        used_decoy_count += 1
        if (used_decoy_count == iDecoyNumber) :
            break
    os.remove(temp_filename)
    output_file = open(output_dir+os.sep+sInternalID_ori+".txt", "w")
    output_file.write(ori_line+"\n")
    decoy_count = 0
    for each_decoy_inchi in decoy_list :
        sInternalID_decoy = "Decoy_"+str(decoy_count)+"_"+sInternalID_ori
        decoy_count += 1
        output_file.write(sInternalID_decoy+"\t"+each_decoy_inchi)
        for i in range(2, len(compound_info_list)):
            output_file.write("\t"+compound_info_list[i])
        output_file.write("\n")
    output_file.close()
    

def GenerateNewDB_mp(input_filename, output_dir, iDecoyNumber) :
    input_file = open(input_filename)
    mypool = Pool(processes=4)
    for each_line in input_file:
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith( "#")) :
            continue
        compound_info_list  = each_line.split("\t")
        sInternalID_ori = compound_info_list[0]
        inchi_ori  = compound_info_list[1]
        result = mypool.apply_async(GenerateDecoy_mp, (sInternalID_ori, inchi_ori, iDecoyNumber, each_line, output_dir, compound_info_list))
    input_file.close()
    mypool.close()
    mypool.join()
    if result.successful():
        print "successful"
        
def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        input_filename ,  output_dir, iDecoyNumber = parse_options(argv)  

    GenerateNewDB_mp(input_filename, output_dir, iDecoyNumber)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()

