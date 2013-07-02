#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def parse_options(argv):

    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                             	     "input-file",
	                			     "output-file"])

    output_filename = ""
    input_filename  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-file, -o output-file"
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-file"):
            output_filename = value

    if (input_filename == "") :
        print "Please specify -i"
        sys.exit(1)
    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_shortern" + inputFileNameExt
    return (input_filename, output_filename)


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

def RemoveTab (sOriginalStr) :
    
    if (sOriginalStr.find("\t") >= 0) :
        print "find tap on", sOriginalStr
        sModifiedStr = sOriginalStr.replace("\t"," ")
    else :
        sModifiedStr = sOriginalStr
    return sModifiedStr

def ExtractFeatureValue (sFeature) :
    sep_pos = sFeature.find(" - ")
    sValue  = sFeature[sep_pos+3:]
    sValue  = sValue.strip()
    sValue  = RemoveTab(sValue)
    return sValue

def CalCompoundMass(sInchi) :
    #print sInchi
    current_mol  = Chem.MolFromInchi(sInchi)
    current_mass = Descriptors.ExactMolWt(current_mol)
    return current_mass

def parseInputFile(input_filename, output_filename) :
    input_file = open(input_filename)
    output_file= open(output_filename, "w")

    sUNIQUE_ID   = ""
    sInchi       = ""
    sDBLink      = "NA"
    sCommonName  = "NA"
    iCompoundNum = 0

    for each_line in input_file :
        each_line = each_line.strip()
        if each_line == "//" :
            iCompoundNum += 1
            sCurrentInternalID = "MetaCyC_"+str(iCompoundNum)
            if ((sUNIQUE_ID == "") or (sInchi == "") or (sDBLink == "")) :
                
                print "one feature of", sUNIQUE_ID, sInchi, "is missed"
            else :
                output_file.write(sCurrentInternalID+"\t")
                output_file.write(sInchi+"\t"+str(dMass)+"\t")
                output_file.write(sCommonName+"\t")
                if (sDBLink != "NA") :
                    output_file.write(sDBLink[3:]+"\n")
                else:
                    output_file.write("NA\n")
             #   output_file.write("//\n")
            sUNIQUE_ID = ""
            sInchi     = ""
            sDBLink    = "NA"
            sCommonName= "NA"
        if each_line.startswith("UNIQUE-ID -") :
            sUNIQUE_ID = ExtractFeatureValue(each_line)
        if each_line.startswith("INCHI -"):
            sInchi     = ExtractFeatureValue(each_line)
            dMass      = CalCompoundMass(sInchi)
        if each_line.startswith("DBLINKS - (PUBCHEM "):
            sDBLink   += "&"+ExtractFeatureValue(each_line)
        if each_line.startswith("COMMON-NAME -") :
            sCommonName = ExtractFeatureValue(each_line)

    input_file.close()
    output_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [input_filename ,  output_filename] = parse_options(argv)  

    parseInputFile(input_filename, output_filename)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()











