#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:c:o:",
                                    ["help",
                                     "working-dir",
				                     "output-dir",])


    # Default working dir and config file
    working_dir = "./"
    output_dir = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory -o outputfile"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-dir"):
            output_dir = value

    SMILES_filename_list = get_file_list_with_ext(working_dir, ".txt")
    
    if (output_dir== "") :
        output_dir = working_dir 

    return [SMILES_filename_list ,  output_dir]

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

def parseMassbankFile(currentInchiFileName, output_file) :
    bPeakBegin = False
    allPeaks_list = []
    current_Inchi_file = open(currentInchiFileName)

    for each_line in current_Inchi_file :    
        output_file.write(each_line)
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

def ExactBondsInfo(FragmentBonds_list) :
    iNumBond = len(FragmentBonds_list)
    if (iNumBond == 0) :
        sBondType = "NA"
    else :
        sBondType = ""
        for current_bond in FragmentBonds_list :
            sBeginAtom = current_bond.GetBeginAtom().GetSymbol()
            sEndAtom   = current_bond.GetEndAtom().GetSymbol()
            sType      = str(current_bond.GetBondType())
            if (sType == "SINGLE"):
                sRepresentType = "-"
            elif (sType == "DOUBLE"):
                sRepresentType = "="
            elif (sType == "TRIPLE"):
                sRepresentType = "~"
            elif (sType == "AROMATIC"):
                sRepresentType = "*"
            else :
                sRepresentType = "("+sType+")"
            sBondType += sBeginAtom + sRepresentType + sEndAtom + ","
        sBondType = sBondType[:-1]
    sAllBond = str(iNumBond)+"\t"+sBondType
    return sAllBond

def NewMatchedFragment (current_peakmatches_list, sCurrent_mz_offset, sCurrent_sFragmentFormula, sCurrent_smiles ) :
    bNewFragment = True
    for each_storedfragment in current_peakmatches_list :
        if (each_storedfragment[0] == sCurrent_mz_offset) :
            if (each_storedfragment[1] == sCurrent_sFragmentFormula) :
                if (each_storedfragment[2] == sCurrent_smiles) :
                    bNewFragment = False
                    break
    return bNewFragment

def MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, FragmentBonds_list) :
    z_list = [1] 
    mz_windows_list = [-1, 0, 1, 2]
    dMass_Tolerance_Fragment_Ions = 0.01
    dHMass = 1.007825
    bFindPeak = False
    bNewFragment = True
    sBondInfo = ExactBondsInfo(FragmentBonds_list)
    for i in range(len(allPeaks_list)) :
        each_peak = allPeaks_list[i]
        dMeasuredMZ = each_peak[0]
        for current_z in z_list :
            for current_mz_offset in mz_windows_list :
                mzdiff =math.fabs ((current_dMass + current_mz_offset*dHMass)/current_z  - dMeasuredMZ)
                if mzdiff <= dMass_Tolerance_Fragment_Ions :
                    dErrorDa = dMeasuredMZ - ((current_dMass+current_mz_offset*dHMass)/current_z)
                    bNewFragment = NewMatchedFragment (peakmatch_list[i], str(current_mz_offset), str(current_sFragmentFormula),str(current_smiles) )
                    if (bNewFragment) :
                        peakmatch_list[i].append([str(current_mz_offset),str(current_sFragmentFormula),str(current_smiles),sBondInfo,dErrorDa])
                    bFindPeak = True
                    break
    return bFindPeak, bNewFragment


def DumpOneFragment(current_fragment_mol, FragmentBonds_list) :
    current_dMass    = Descriptors.ExactMolWt(current_fragment_mol)
    current_sFragmentFormula = AllChem.CalcMolFormula(current_fragment_mol)
    sBondsTypes = "{"
    for current_bond in FragmentBonds_list :
        if (sBondsTypes != "{") :
            sBondsTypes += ","
        sBondsTypes += str(current_bond.GetBondType())
    sBondsTypes += "}"
    #current_SanitizedMol = Chem.SanitizeMol(current_fragment_mol)
    #current_inchi =  Chem.MolToInchi(current_fragment_mol)
    current_smiles=  Chem.MolToSmiles (current_fragment_mol)
    #print current_dMass, current_sFragmentFormula, sBondsTypes, current_smiles
    return current_dMass, current_sFragmentFormula, current_smiles

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

def RemoveBonds(current_mol, bonds_list) :
   # print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    em = Chem.EditableMol(current_mol)
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
    #print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    return current_fragments_list, bValidOperation


def processKid(current_editable_mol, current_removebond_list, iCurrent_depth, allPeaks_list, peakmatch_list) :
    #print current_editable_mol
    current_mol = current_editable_mol.GetMol()
    current_dMass, current_sFragmentFormula, current_smiles=DumpOneFragment(current_mol, current_removebond_list)
    bFindPeak, bNewFragment=MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, current_removebond_list)
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(current_mol)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []

    new_item = [[current_editable_mol, current_removebond_list, iCurrent_depth], current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid]
    return new_item


def TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, iDepth, peakmatch_list) :
    root_node = [Chem.EditableMol(current_mol),[],0] # editable_mol,list of list of removed bonds, depth
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(root_node[0].GetMol())
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []
    storedNodes = [[root_node, current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid]]
    while (len(storedNodes) > 0) :
        if (len(storedNodes[-1][3]) > 0) : # unprocessed kid
            new_item = processKid(storedNodes[-1][3][0][0], storedNodes[-1][3][0][1], storedNodes[-1][0][2]+1, allPeaks_list, peakmatch_list) 
            del storedNodes[-1][3][0]
            if (new_item[0][2] < iDepth) :
                storedNodes.append(new_item)
        elif (len(storedNodes[-1][1]) > 0) : # linear bond
            remove_bond = storedNodes[-1][1][0]
            current_fragments_list, bValidOperation = RemoveBonds(storedNodes[-1][0][0].GetMol(), [remove_bond]) 
            if (bValidOperation) :
                FragmentBonds_list  =  list( storedNodes[-1][0][1] )
                FragmentBonds_list.append(remove_bond)
                for i in range(2) :
                    storedNodes[-1][3].append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list])
            else :
                print "wrong!"
                sys.exit(1)
            del storedNodes[-1][1][0]
        elif (len(storedNodes[-1][2]) > 0) : # ring bonds
            remove_first_bond  = storedNodes[-1][2][0][0]
            remove_second_bond = storedNodes[-1][2][0][1]
            current_fragments_list, bValidOperation = RemoveBonds(storedNodes[-1][0][0].GetMol(), [remove_first_bond, remove_second_bond])
            if (bValidOperation) :
                FragmentBonds_list  =  list( storedNodes[-1][0][1] )
                FragmentBonds_list.append(remove_first_bond)
                FragmentBonds_list.append(remove_second_bond)
                for i in range(2) :
                    storedNodes[-1][3].append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list])
            del storedNodes[-1][2][0]
        else :
            del storedNodes[-1]


def ExhaustBonds(sInchiInfo, allPeaks_list, peakmatch_list) :
    current_mol = Chem.MolFromInchi(sInchiInfo)
    current_sFragmentFormula = AllChem.CalcMolFormula(current_mol)
    current_smiles = Chem.MolToSmiles(current_mol)
    print Descriptors.ExactMolWt(current_mol), "NULL", AllChem.CalcMolFormula(current_mol)
    MapMass(Descriptors.ExactMolWt(current_mol), allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, [])
    iBondsNum = current_mol.GetNumBonds()
    #TreeLikeBreakBonds(current_mol, iBondsNum, allPeaks_list, 3, peakmatch_list)
    TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, 3, peakmatch_list)

def HandleInchi(currentInchiFileName, output_file) :
    #output_file.write(">\t"+currentInchiFileName+"\n")
    sInchiInfo, allPeaks_list = parseMassbankFile(currentInchiFileName, output_file)
    output_file.write("\n")
    output_file.write("peak\tm/z\tint\trel_int\tH_added\tformula\tSMILES\tnum_bonds_cleaved\ttypes_bonds_cleaved\terror_Da\n")
    peakmatch_list = [[] for each_peak in allPeaks_list]
    ExhaustBonds(sInchiInfo, allPeaks_list, peakmatch_list)
    for i in range(len(allPeaks_list)):
        current_peak = allPeaks_list[i]
        if (len(peakmatch_list[i]) == 0 ) :
            output_file.write(str(i+1))
            for peak_item in current_peak :
                output_file.write("\t" + str(peak_item))
            output_file.write("\tNA\tNA\tNA\tNA\tNA\tNA\n")
        else :
            for each_peakmatch in peakmatch_list[i] :
                output_file.write(str(i+1))
                for peak_item in current_peak :
                    output_file.write("\t" + str(peak_item))
                for each_data in each_peakmatch :
                     output_file.write("\t" + str(each_data))
                output_file.write("\n")   


def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [Inchi_filename_list ,  output_dir] = parse_options(argv)  
    #output_file = open(outputFileName, "w")
    for eachInchiFileName in Inchi_filename_list :
        outputFileName = output_dir+os.path.basename(eachInchiFileName)+"output.txt"
        output_file = open(outputFileName, "w")
        print eachInchiFileName
        HandleInchi(eachInchiFileName, output_file)
        print "***********************************************************"
        #os.remove(eachInchiFileName)
        output_file.close()



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()