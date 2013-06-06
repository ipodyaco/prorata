#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

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
    bPeakBegin = False
    allPeaks_list = []
    current_Inchi_file = open(currentInchiFileName)

    for each_line in current_Inchi_file :    
        each_line = each_line.strip()
        if each_line.startswith("CH$IUPAC:") :
            s_chemical_structure = each_line.split(" ")[1]
            #print s_chemical_structure
        if each_line.startswith("PK$NUM_PEAK:") :
            iPeakNum =int (each_line.split(":")[1])
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

    current_Inchi_file.close()
    #print allPeaks_list 
    return s_chemical_structure, allPeaks_list

def MapMass(current_dMass, allPeaks_list) :
    z_list = [1] 
    mz_windows_list = [-2, -1, 0, 1, 2]
    dMass_Tolerance_Fragment_Ions = 0.01


def BreakOneBond(current_mol, iBondsNum, allPeaks_list) :
    for i in range(iBondsNum) :
        current_editable_mol = Chem.EditableMol(current_mol)
        current_bond  = current_mol.GetBondWithIdx(i)
        idx_beginAtom = current_bond.GetBeginAtomIdx()
        idx_endAtom   = current_bond.GetEndAtomIdx()
        current_editable_mol.RemoveBond(idx_beginAtom, idx_endAtom)
        current_modified_mol = current_editable_mol.GetMol()
        current_fragments    = Chem.GetMolFrags(current_modified_mol, asMols=True, sanitizeFrags=False)
        #current_fragments    = Chem.GetMolFrags(current_modified_mol)
        if (len(current_fragments) > 1) :
            for each_fragment in current_fragments :
                #current_sInchi = Chem.MolToInchi(each_fragment)
                current_dFragmentFormula = AllChem.CalcMolFormula(each_fragment)
                current_dMass  = Descriptors.ExactMolWt(each_fragment)
                current_sBondType = current_bond.GetBondType()
                print current_dMass, current_sBondType, current_dFragmentFormula
        #        print each_fragment
                MapMass(current_dMass, allPeaks_list)

def RemoveOneBond(original_mol, current_editable_mol, bond_idx) :
    current_bond  = original_mol.GetBondWithIdx(bond_idx)
    idx_beginAtom = current_bond.GetBeginAtomIdx()
    idx_endAtom   = current_bond.GetEndAtomIdx()
    current_editable_mol.RemoveBond(idx_beginAtom, idx_endAtom)
    current_modified_mol = current_editable_mol.GetMol()
    current_fragments_list = Chem.GetMolFrags(current_modified_mol, asMols=True, sanitizeFrags=False)
    return current_editable_mol, current_fragments_list

def DumpFragments(current_mol, fragment_list) :
    for each_fragment_info in fragment_list :
        current_fragment = each_fragment_info[0]
        current_bond_idx = each_fragment_info[1]
        current_bond     = current_mol.GetBondWithIdx(current_bond_idx)
        current_sBondType = current_bond.GetBondType()
        current_dFragmentFormula = AllChem.CalcMolFormula(current_fragment)
        current_dMass    = Descriptors.ExactMolWt(current_fragment)
        print current_dMass, current_sBondType, current_dFragmentFormula


def BreakBonds(current_mol, iBondsNum, allPeaks_list, iMaxBreakBondNum) :
    seed_mol_list = []
    new_seed_mol_list = []
    fragment_list = []
    iRealMaxBreakBondNum = min(iBondsNum, iMaxBreakBondNum)
    for iBreakBondNum in range(1, iRealMaxBreakBondNum+1) :
        if (iBreakBondNum == 1) :
            for BreakBond_idx in range(iBondsNum) :
                #bGenerateNewFragment = False
                current_editable_mol = Chem.EditableMol(current_mol)
                current_updated_editable_mol, current_fragments_list = RemoveOneBond(current_mol, current_editable_mol, BreakBond_idx)
                seed_mol_list.append([current_updated_editable_mol, BreakBond_idx])
                if (len(current_fragments_list) > 1) : 
                    for each_fragment in current_fragments_list :
                        fragment_list.append([each_fragment, BreakBond_idx])
        else :
            new_seed_mol_list = []
            for each_seed in seed_mol_list :
                seed_mol = each_seed[0].GetMol()
                seed_fragment_num = len(Chem.GetMolFrags(seed_mol, asMols=True, sanitizeFrags=False))
                for BreakBond_idx in range(each_seed[1]+1, iBondsNum) :
                    current_editable_mol = each_seed[0]
                    current_updated_editable_mol, current_fragments_list = RemoveOneBond(current_mol, current_editable_mol, BreakBond_idx)
                    new_seed_mol_list.append([current_updated_editable_mol, BreakBond_idx])
                    if (len(current_fragments_list) > seed_fragment_num) :
                        for each_fragment in current_fragments_list :
                            fragment_list.append([each_fragment, BreakBond_idx])
            seed_mol_list = new_seed_mol_list
            new_seed_mol_list =[]
    DumpFragments(current_mol, fragment_list) 

def ExhaustBonds(sInchiInfo, allPeaks_list) :
    current_mol = Chem.MolFromInchi(sInchiInfo)
    print Descriptors.ExactMolWt(current_mol), "NULL", AllChem.CalcMolFormula(current_mol)
    iBondsNum = current_mol.GetNumBonds()
    #BreakOneBond(current_mol, iBondsNum, allPeaks_list)
    BreakBonds(current_mol, iBondsNum, allPeaks_list, 2)
    

def HandleInchi(currentInchiFileName, output_file) :
    output_file.write(">\t"+currentInchiFileName+"\n")
    sInchiInfo, allPeaks_list = parseMassbankFile(currentInchiFileName)
    ExhaustBonds(sInchiInfo, allPeaks_list)
    


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
