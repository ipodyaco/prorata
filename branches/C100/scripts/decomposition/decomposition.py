#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
import itertools

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


def DumpOneFragment(current_fragment_mol, FragmentBonds_list) :
    current_dMass    = Descriptors.ExactMolWt(current_fragment_mol)
    current_sFragmentFormula = AllChem.CalcMolFormula(current_fragment_mol)
    sBondsTypes = "{"
    for current_bond in FragmentBonds_list :
        if (sBondsTypes != "{") :
            sBondsTypes += ","
        sBondsTypes += str(current_bond.GetBondType())
    sBondsTypes += "}"
    current_inchi =  Chem.MolToInchi(current_fragment_mol)
    print current_dMass, current_sFragmentFormula, sBondsTypes, current_inchi

def ClassifyBonds(current_mol) :
    ring_bonds_list   = []
    linear_bonds_list = []
    iBondsNum = current_mol.GetNumBonds() 
    for i in range(iBondsNum) :
        current_bond = current_mol.GetBondWithIdx(i)
        if current_bond.IsInRing()  :
            ring_bonds_list.append(current_bond)
        else :
            linear_bonds_list.append(current_bond)
    return ring_bonds_list, linear_bonds_list

def RemoveBonds(current_editable_mol, bonds_list) :
    em = current_editable_mol
    for each_removable_bond in bonds_list :
        idx_beginAtom = each_removable_bond.GetBeginAtomIdx()
        idx_endAtom   = each_removable_bond.GetEndAtomIdx()
        em.RemoveBond(idx_beginAtom, idx_endAtom)
    current_modified_mol = em.GetMol()
    current_fragments_list = Chem.GetMolFrags(current_modified_mol, asMols=True, sanitizeFrags=False)
    if len(current_fragments_list) == 2 :
        bValidOperation = True
    else :
        bValidOperation = False
    return current_fragments_list, bValidOperation


def TreeLikeBreakBonds(current_mol, iBondsNum, allPeaks_list, iDepth) :
    FragmentTree_list = [[Chem.EditableMol(current_mol), [], [], -1]] # editable_mol, list of list of removed bonds, children, father Idx 
    iCurrentGrandpaIdx= -1
    iCurrentFatherIdx = 0
    iCurrentChindIdx  = 0
    for iCurrentDepth in range (iDepth) :
        if (iCurrentGrandpaIdx < 0) : # at root
            iCurrent_FatherIdx_list = [0]
        else :
            iCurrent_FatherIdx_list = FragmentTree_list[iCurrentGrandpaIdx][2]
        for iCurrentFatherIdx in iCurrent_FatherIdx_list :
            FatherFragmentInfo_list = FragmentTree_list [iCurrentFatherIdx]
            current_editable_mol    = FatherFragmentInfo_list[0]
            #Classify Bonds
            current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(FatherFragmentInfo_list[0].GetMol())
            
            for each_bond in current_linear_bonds_list :
                current_fragments_list, bValidOperation = RemoveBonds(current_editable_mol, [each_bond])
                if (bValidOperation) :
                    for i in range(2) :
                        iCurrentFragmentNum =  len(FragmentTree_list)
                        FragmentBonds_list  =  list( FragmentTree_list [iCurrentFatherIdx][1] )
                        print FragmentBonds_list
                        FragmentTree_list [iCurrentFatherIdx][2].append(iCurrentFragmentNum)
                        FragmentBonds_list.append(each_bond)
                        FragmentTree_list.append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list, [], iCurrentFatherIdx])
                        DumpOneFragment(current_fragments_list[i], FragmentBonds_list)
                else :
                    print "wrong!"
                    sys.exit(1)

            removable_bonds_iter = itertools.combinations(current_ring_bonds_list, 2)
            for first_bond, second_bond in removable_bonds_iter :
                current_fragments_list, bValidOperation = RemoveBonds(current_editable_mol, [first_bond, second_bond])
                if (bValidOperation) :
                    for i in range(2) :
                        iCurrentFragmentNum =  len(FragmentTree_list) 
                        FragmentBonds_list  =  list( FragmentTree_list [iCurrentFatherIdx][1] )
                        FragmentTree_list [iCurrentFatherIdx][2].append(iCurrentFragmentNum)
                        FragmentBonds_list.append(first_bond)
                        FragmentBonds_list.append(second_bond)
                        FragmentTree_list.append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list, [], iCurrentFatherIdx])
                        DumpOneFragment(current_fragments_list[i], FragmentBonds_list)


def ExhaustBonds(sInchiInfo, allPeaks_list) :
    current_mol = Chem.MolFromInchi(sInchiInfo)
    print Descriptors.ExactMolWt(current_mol), "NULL", AllChem.CalcMolFormula(current_mol)
    iBondsNum = current_mol.GetNumBonds()
    TreeLikeBreakBonds(current_mol, iBondsNum, allPeaks_list, 3)
    

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
        print "***********************************************************"

    output_file.close()



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
