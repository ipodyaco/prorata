#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hw:o:m:",
                                    ["help",
                                     "working-dir",
				                     "output-filename",
                                     "mode"])


    # Default working dir and config file
    working_dir       = ""
    output_filename   = ""
    mode              = 0

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w working-dir -o output_filename -m mode(1 for [M+H]+ and -1 for [M-H]-)"
            sys.exit(0)
        if option in ("-w", "--working-dir"):
            working_dir  = value
            if (working_dir[-1] != os.sep) :
                working_dir += os.sep
        if option in ("-o", "--output-filename"):
            output_filename   = value
        if option in ("-m", "--mode") :
            mode = int(value)
    if ( working_dir == ""  ) :
        print "please specify working dir"
        sys.exit(1)
    if ((mode != 1) and (mode != -1)) :
        print "wrong mode"
        sys.exit(1)

    if (output_filename == "") :
        output_filename = working_dir+"PubChem_Massbank.txt"
    
    input_filename_list = get_file_list_with_ext(working_dir, "txt")

    return [input_filename_list, output_filename, working_dir, mode]

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

def VerifyStructure(s_chemical_structure) :
    #print s_chemical_structure
    current_mol = Chem.MolFromInchi(s_chemical_structure)
    dCalculatedMass = Descriptors.ExactMolWt(current_mol)
    if (current_mol == None) :
        bLegalStructure = False
    else :
        bLegalStructure = True
    return bLegalStructure, dCalculatedMass

def ReadRealHit(realhit_filename):
    
    dProtonMass = 1.007825 # proton mass
    bPeakBegin = False
    allPeaks_list = []
    precursor_mz = -1
    realhit_file = open(realhit_filename)
    sPubChemId   = ""
    sMSType      = ""
    sPrecursorType = ""
    sCompoundName  = ""

    for each_line in realhit_file :
        each_line = each_line.strip()
        if each_line.startswith("CH$IUPAC:") :
            s_chemical_structure = each_line.split(" ")[1].strip()
        if each_line.startswith("PK$NUM_PEAK:") :
            iPeakNum =int (each_line.split(":")[1])
        if each_line.startswith("MS$FOCUSED_ION: PRECURSOR_M/Z ") :
            precursor_mz = float(each_line.split("PRECURSOR_M/Z")[1])
        if each_line.startswith("CH$LINK: PUBCHEM CID:") :
            sPubChemId = each_line.split("PUBCHEM CID:")[1].strip()
        if each_line.startswith("ACCESSION:") :
            sAccession = each_line.split(":")[1].strip()
        if each_line.startswith("AC$MASS_SPECTROMETRY: MS_TYPE "):
            sMSType = each_line.split("MS_TYPE")[1].strip()
        if each_line.startswith("MS$FOCUSED_ION: PRECURSOR_TYPE ") :
            sPrecursorType = each_line.split("PRECURSOR_TYPE")[1].strip()
        if each_line.startswith("CH$NAME:") :
            sCompoundName += "&"+each_line.split(":")[1].strip()
        if each_line.startswith("PK$PEAK:") :
            bPeakBegin = True
            continue
        if bPeakBegin :
            sPeakInfo_list = each_line.split(" ")
            currentPeak = []
            # Peak format m/z int rel.int
            for i in range(3) :
                if (i<2) :
                    currentPeak.append(sPeakInfo_list[i].strip())
                else :
                    currentPeak.append(sPeakInfo_list[i].strip())
            allPeaks_list.append(currentPeak)
            iPeakNum -= 1
            if (iPeakNum == 0) :
                bPeakBegin = False

    realhit_file.close()
    current_mode = 0
    bLegalStructure,dCalculatedMass = VerifyStructure(s_chemical_structure)
    if (sPrecursorType == "[M+H]+") :
        current_mode = 1
        massdiff = math.fabs(precursor_mz-dProtonMass - dCalculatedMass)
        if ((massdiff >= 0.02) and (sMSType == "MS2")):
            print realhit_filename, massdiff, "1"
    if (sPrecursorType == "[M-H]-") :
        current_mode = -1
        massdiff = math.fabs(precursor_mz+dProtonMass - dCalculatedMass)
        if ((massdiff >= 0.02) and (sMSType == "MS2")):
            print realhit_filename, massdiff, "-1"
    return precursor_mz - current_mode*dProtonMass, allPeaks_list, s_chemical_structure, sPubChemId, sAccession, sMSType, bLegalStructure, sPrecursorType, current_mode, sCompoundName

def extractMassBankFile(input_filename, output_file, mode) :

    dPrecursorEM, sPeak_list, sInichi, sPubChemId, sAccession, sMSType, bLegalStructure, sPrecursorType, current_mode, sCompoundName = ReadRealHit(input_filename)

#    if "CE000438" in  input_filename:
#        print bLegalStructure

    if (sPubChemId == "") :
        sPubChemId = "NA"
    if (sCompoundName == "") :
        sCompoundName = "NA"
    if ((sMSType == "MS2") and (bLegalStructure)and (mode*current_mode == 1)):
        if (dPrecursorEM < 0) :
            print input_filename, "no precursor"
        output_file.write("*\t"+sAccession+"\t"+sPubChemId+"\t"+str(dPrecursorEM)+"\t"+sInichi+"\t"+sCompoundName+"\n")
        for each_peak in sPeak_list :
            output_file.write("+"+each_peak[0]+" "+each_peak[1]+"\n")
        output_file.write("//\n")
    elif not (bLegalStructure) :
        print input_filename, "illegal structure", sInichi
#    else :
#        print input_filename, sMSType

    return sPrecursorType

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [input_filename_list, output_filename, working_dir, mode] = parse_options(argv)  
        #badcase_dir = working_dir+"nopubchemid"
        #os.system("mkdir " + badcase_dir)
        output_file = open(output_filename, "w")
        
        sPrecursorType_list = []
        iPrecursorType_list =[]

        for each_input_filename in input_filename_list :
            sPrecursorType = extractMassBankFile(each_input_filename, output_file, mode)
            if (sPrecursorType in sPrecursorType_list) :
                iPrecursorType_list[sPrecursorType_list.index(sPrecursorType)] += 1
            else :
                sPrecursorType_list.append(sPrecursorType)
                iPrecursorType_list.append(1)
        for i in range(len(iPrecursorType_list)) :
            print sPrecursorType_list[i], iPrecursorType_list[i]
        output_file.close()

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



