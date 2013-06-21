#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hm:p:o:e:",
                                    ["help",
                                     "realhit-filename",
                                     "compound-filename",
				                     "output-filename",
                                     "energy-filename"])


    # Default working dir and config file
    realhit_filename  = ""
    compound_filename = ""
    output_filename   = ""
    energy_filename   = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-m realhit_filename -p compound_filename -o output_filename -e energy_filename"
            sys.exit(0)
        if option in ("-m", "--realhit-filename"):
            realhit_filename  = value
        if option in ("-p", "--compound-filename"):
            compound_filename = value
        if option in ("-o", "--output-filename"):
            output_filename   = value
        if option in ("-e", "--energy-filename"):
            energy_filename   = value

    if ((realhit_filename == "") or (compound_filename == "") or (output_filename == "") or (energy_filename == "")) :
        print "please specify realhit file, compound file, output file, and energy file"
        sys.exit(1)
    
    return [realhit_filename, compound_filename, output_filename, energy_filename]

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


def ReadEnergyFile(energy_filename) :
    
    Bond_Symbols_list = ["-","=","~"] 
    sEnergy_Bond_dict = {}
    energy_file = open(energy_filename)

    for each_line in energy_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line == "//"):
            break
        energy_info_list = each_line.split(" ")
        sCurrent_Bond    = energy_info_list[0]
        dCurrent_Energy  = float(energy_info_list[-1])
        sEnergy_Bond_dict[sCurrent_Bond] = dCurrent_Energy
        for sCurrentBondSymbol in Bond_Symbols_list :
            if sCurrentBondSymbol in sCurrent_Bond :
                sAtom_Symol_list = sCurrent_Bond.split(sCurrentBondSymbol)
                dMirror_Bond     = sAtom_Symol_list[-1] + sCurrentBondSymbol + sAtom_Symol_list[0]
                sEnergy_Bond_dict[dMirror_Bond] = dCurrent_Energy
                continue

    energy_file.close()
    return sEnergy_Bond_dict

def ReadCompoundFile(compound_filename) :
    compound_list  = [] # each entry: ID, Inchi, Mass, and DBLink
    compound_file = open(compound_filename)
    for each_line in compound_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        current_compound_info = each_line.split("\t")
        if (len(current_compound_info) != 4) :
            print "illegal compound", each_line
            sys.exit(1)
        else:
            compound_list.append([current_compound_info[0], current_compound_info[1], float(current_compound_info[2]), current_compound_info[3]])

    compound_file.close()

    compound_list.sort(key=lambda e:e[2]) # sort based on mass
    return compound_list

def ReadRealHit(realhit_filename):

    bPeakBegin = False
    allPeaks_list = []
    precursor_mz = -1
    realhit_file = open(realhit_filename)
    
    for each_line in realhit_file :
        each_line = each_line.strip()
        if each_line.startswith("CH$IUPAC:") :
            s_chemical_structure = each_line.split(" ")[1]
        if each_line.startswith("PK$NUM_PEAK:") :
            iPeakNum =int (each_line.split(":")[1])
        if each_line.startswith("MS$FOCUSED_ION: PRECURSOR_M/Z ") :
            precursor_mz = float(each_line.split("PRECURSOR_M/Z")[1])
        if each_line.startswith("PK$PEAK:") :
            bPeakBegin = True
            continue
        if bPeakBegin :
            sPeakInfo_list = each_line.split(" ")
            currentPeak = []
            # Peak format m/z int rel.int
            for i in range(3) :
                if (i<2) :
                    currentPeak.append(float(sPeakInfo_list[i]))
                else :
                    currentPeak.append(int(sPeakInfo_list[i]))
            allPeaks_list.append(currentPeak)
            iPeakNum -= 1
            if (iPeakNum == 0) :
                bPeakBegin = False

    realhit_file.close()
    return precursor_mz, allPeaks_list

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [realhit_filename, compound_filename, output_filename, energy_filename] = parse_options(argv)  

    sEnergy_Bond_dict = ReadEnergyFile(energy_filename)
    #print sEnergy_Bond_dict.items()
    Compound_list    = ReadCompoundFile(compound_filename)
    #print Compound_list[30]
    precursor_mz, allPeaks_list = ReadRealHit(realhit_filename)
    #print precursor_mz, allPeaks_list     

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()


