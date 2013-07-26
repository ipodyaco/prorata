#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

import weightedscore


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
        if ((each_line == "") or (each_line.startswith("#"))) :
            continue
        current_compound_info = each_line.split("\t")
        if (len(current_compound_info) != 5) :
            print "illegal compound", each_line
            sys.exit(1)
        else:
            compound_list.append([current_compound_info[0], current_compound_info[1], float(current_compound_info[2]), current_compound_info[3]])

    compound_file.close()

    compound_list.sort(key=lambda e:e[2]) # sort based on mass
    return compound_list

def ReadRealHit(realhit_filename):
    
    dProtonMass = 1.007825 # proton mass
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
    return precursor_mz-dProtonMass, allPeaks_list, s_chemical_structure

def BinarySearch_Upper(target_list, bounder_value):
    ele_num = len(target_list)
    upper_index = ele_num-1
    if (target_list[upper_index] < bounder_value):
        return upper_index
    lower_index = 0
    if (target_list[lower_index] > bounder_value) :
        return -1
    while ((upper_index - lower_index)>1):
        middle_index = (lower_index + upper_index)/2
        if (target_list[middle_index]> bounder_value ) :
            upper_index = middle_index
        else :
            lower_index = middle_index

    return lower_index
    
def BinarySearch_Lower(target_list, bounder_value):
    ele_num = len(target_list)
    upper_index = ele_num-1
    if (target_list[upper_index] < bounder_value):
        return -1
    lower_index = 0
    if (target_list[lower_index] > bounder_value) :
        return lower_index
    while ((upper_index - lower_index)>1):
        middle_index = (lower_index + upper_index)/2
        if (target_list[middle_index]>= bounder_value ) : # !!!!
            upper_index = middle_index
        else :
            lower_index = middle_index
    return upper_index

def GetRelatedCompound(Compound_list, precursor_mz, precursor_accuracy):
    compound_mass_list = [each_compound[2] for each_compound in Compound_list]
    upper_compound_mass = precursor_mz + precursor_accuracy
    lower_compound_mass = precursor_mz - precursor_accuracy
    upper_compound_index=BinarySearch_Upper(compound_mass_list, upper_compound_mass)
    #print upper_compound_index, compound_mass_list[upper_compound_index], upper_compound_mass, compound_mass_list[upper_compound_index+1]
    lower_compound_index=BinarySearch_Lower(compound_mass_list, lower_compound_mass)
    #print lower_compound_index, compound_mass_list[lower_compound_index-1], lower_compound_mass, compound_mass_list[lower_compound_index] 
    QueryCompound_list = []
    if ((lower_compound_index != -1) and (upper_compound_index != -1)) :
        QueryCompound_list = Compound_list[lower_compound_index : upper_compound_index+1]
    return QueryCompound_list


def NormalizeIntensity(allPeaks_list) :
    max_intensity = max([each_peak[1] for each_peak in allPeaks_list])
    for i in range(len(allPeaks_list)) :
        allPeaks_list[i][1] = allPeaks_list[i][1]/max_intensity
    return allPeaks_list

def RankScores(Compound_Scores_list, output_filename, bSpectrumDetails, sRealInchi) :
    sCoreName = os.path.splitext(os.path.basename(output_filename))[0]
    if (sCoreName.endswith("output")) :
        sCoreName = os.path.splitext(sCoreName)[0]
    Compound_Scores_list.sort(key=lambda e:e[0], reverse=True)
    ID_list = [each_compound_score[3] for each_compound_score in Compound_Scores_list ]
    realhit_index = ID_list.index(sRealInchi) 
    realhit_score = Compound_Scores_list[realhit_index][0]
   # print realhit_score
    output_str = "" 
    iRealhitRank = 0
    current_rank = 1
    for i in range(len(Compound_Scores_list)):
        dCurrentScore = Compound_Scores_list[i][0]
        if (i==0) :
            current_rank = 1
            dPreviousScore = dCurrentScore
        else :
            if (dCurrentScore < dPreviousScore):
                current_rank = i+1
            dCurrentScore = dPreviousScore
        if (i == realhit_index) :
            iRealhitRank = current_rank
            #print iRealhitRank
        output_str +=str(current_rank)
        for each_item in Compound_Scores_list[i][:-1] :
            output_str += "\t"+str(each_item)
        output_str += "\n"
        if (bSpectrumDetails) :
            for each_annontation in Compound_Scores_list[i][-1] :
                output_str += each_annontation + "\n"
   # print sCoreName
   # print ">"+str(iRealhitRank)+"\t"+sCoreName+"\n"
    output_file = open(output_filename, "w")
    output_file.write(">"+str(iRealhitRank)+"\t"+sCoreName+"\n")
    output_file.write(output_str)
    output_file.close()



def score_main(compound_filename, output_filename, energy_filename, bSpectrumDetails, bBreakRing, precursor_type, bRankSum, precursor_mz, allPeaks_list, s_chemical_structure):

    precursor_accuracy = 0.5
    # [realhit_filename, compound_filename, output_filename, energy_filename, bSpectrumDetails, bBreakRing, precursor_type] = parse_options(argv)  
#    print output_filename    
    Compound_Scores_list = []
    sEnergy_Bond_dict = ReadEnergyFile(energy_filename)
    Compound_list    = ReadCompoundFile(compound_filename)
    #precursor_mz, allPeaks_list, s_chemical_structure = ReadRealHit(realhit_filename)
    allPeaks_list = NormalizeIntensity(allPeaks_list)
    QueryCompound_list = GetRelatedCompound(Compound_list, precursor_mz, precursor_accuracy)
#    print len(QueryCompound_list)
    for each_compound in QueryCompound_list :
        current_mol = Chem.MolFromInchi(each_compound[1])
        current_fragments_list = Chem.GetMolFrags(current_mol, asMols=True, sanitizeFrags=False)
        if (len(current_fragments_list) != 1) :
            print "filter", each_compound[1]
            continue
        dCurrentScore, dCurrentEnergy, iIdentifiedPeak, sAnnotation_list, sOtherInfo= weightedscore.OwnScore(sEnergy_Bond_dict, allPeaks_list, current_mol, bBreakRing, precursor_type, bRankSum)
#        print dCurrentScore
        Compound_Scores_list.append([dCurrentScore, dCurrentEnergy, each_compound[0], each_compound[1], each_compound[3], iIdentifiedPeak, sOtherInfo,  sAnnotation_list])
#    print len(Compound_Scores_list), output_filename
    RankScores(Compound_Scores_list, output_filename, bSpectrumDetails, s_chemical_structure)

## If this program runs as standalone, then go to main.
#if __name__ == "__main__":
#    main()


