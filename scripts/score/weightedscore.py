#!/usr/bin/python


import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


#global dMass_Tolerance_Fragment_Ions
#dMass_Tolerance_Fragment_Ions = 0.01


def OwnScore(sEnergy_Bond_dict, allPeaks_list, current_mol, bBreakRing) :
#    global dMass_Tolerance_Fragment_Ions
#    dMass_Tolerance_Fragment_Ions = 0.01
    sOtherInfo = "otherinfo"
    dCurrentScore = 0
    dCurrentEnergy = 0
    iIdentifiedPeak= 0
    peakmatch_list = [[] for each_peak in allPeaks_list]
    ExhaustBonds(current_mol, allPeaks_list, peakmatch_list, sEnergy_Bond_dict, bBreakRing)
   # print peakmatch_list
    sAnnotation_list = []
    sAnnotation_list.append("m/z\tnormalized_int\trel_int\tH_added\tformula\tSMILES\types_bonds_cleaved\terror_Da")
    for i in range(len(allPeaks_list)) :
        sCurrentAnnotation   = str(allPeaks_list[i][0])+"\t"+str(allPeaks_list[i][1])+"\t"+str(allPeaks_list[i][2])
        if peakmatch_list[i] :
            iIdentifiedPeak += 1 
            dCurrentEnergy  += peakmatch_list[i][0][5]
            dErrorDa         = peakmatch_list[i][0][4]
            dSubScore        = peakmatch_list[i][0][6]
            dCurrentScore   += dSubScore
            sCurrentAnnotation += "\t" + str(peakmatch_list[i][0][0]) + "\t" + str(peakmatch_list[i][0][1])+"\t"+str(peakmatch_list[i][0][2])
            sCurrentAnnotation += "\t" + str(peakmatch_list[i][0][3]) + "\t" + str(peakmatch_list[i][0][4])
        else :
            sCurrentAnnotation += "\tNA"
        sAnnotation_list.append(sCurrentAnnotation)
    return dCurrentScore, dCurrentEnergy, iIdentifiedPeak, sAnnotation_list, sOtherInfo


def TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, iDepth, peakmatch_list, sEnergy_Bond_dict, energy_denominator, bBreakRing) :
    root_node = [Chem.EditableMol(current_mol),[],0] # editable_mol,list of list of removed bonds, depth
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(root_node[0].GetMol(), bBreakRing)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []
    storedNodes = [[root_node, current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid]]
    while (len(storedNodes) > 0) :
        if (len(storedNodes[-1][3]) > 0) : # unprocessed kid
            new_item = processKid(storedNodes[-1][3][0][0], storedNodes[-1][3][0][1], storedNodes[-1][0][2]+1, allPeaks_list, peakmatch_list, sEnergy_Bond_dict, energy_denominator, bBreakRing) 
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

def ExactBondsInfo(FragmentBonds_list, sEnergy_Bond_dict, energy_denominator) :
    dF       = 1.0
    dBDE     = 0
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
            #elif (sType == "AROMATIC"):
            #    sRepresentType = "*"
            else :
                sRepresentType = "("+sType+")"
            sCurrentBondType = sBeginAtom + sRepresentType + sEndAtom 
            if sCurrentBondType in sEnergy_Bond_dict :
                dCurrentBDE = sEnergy_Bond_dict[sCurrentBondType]
                dBDE += dCurrentBDE
            else :
                dCurrentBDE = 348  #C-C
                dBDE += dCurrentBDE
            sBondType += sCurrentBondType + ","
            dF = dF * ((1/dCurrentBDE)/energy_denominator)
        sBondType = sBondType[:-1]
    sAllBond = str(iNumBond)+"\t"+sBondType
    return sAllBond, dBDE, dF

def mypnorm (dMean, dStandardDeviation, dRandomVariable) :
    dZScore = ( dRandomVariable - dMean ) / dStandardDeviation
    dProbability = 0.5 * math.erfc( -1 * dZScore / math.sqrt( 2.0 ) )
    return dProbability

def SubScore(dIntensity, dErrorDa, dF, dMass_Tolerance_Fragment_Ions, current_mz_offset):
    dSubScore = 0
    dErrorScore = ( 1.0 - mypnorm( 0, ( dMass_Tolerance_Fragment_Ions / 2.0), math.fabs( dErrorDa  ) ) ) * 2.0
    dSubScore = dIntensity * dErrorScore * dF
    return dSubScore

def MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, FragmentBonds_list, sEnergy_Bond_dict, energy_denominator) :
    z_list = [1] 
    mz_windows_list = [-1, 0, 1, 2]
    dMass_Tolerance_Fragment_Ions = 0.01
    dHMass = 1.007825
    bFindPeak = False
    sBondInfo, dBDE, dF = ExactBondsInfo(FragmentBonds_list, sEnergy_Bond_dict, energy_denominator)
    for i in range(len(allPeaks_list)) :
        each_peak = allPeaks_list[i]
        dMeasuredMZ = each_peak[0]
        for current_z in z_list :
            for current_mz_offset in mz_windows_list :
                mzdiff =math.fabs ((current_dMass + current_mz_offset*dHMass)/current_z  - dMeasuredMZ)
                if mzdiff <= dMass_Tolerance_Fragment_Ions :
                    dErrorDa = dMeasuredMZ - ((current_dMass+current_mz_offset*dHMass)/current_z)
                    dSubScore= SubScore(allPeaks_list[i][1], dErrorDa, dF, dMass_Tolerance_Fragment_Ions, current_mz_offset)
                    current_newmatch = [str(current_mz_offset),str(current_sFragmentFormula),str(current_smiles),sBondInfo,dErrorDa,dBDE, dSubScore]
                    if (not peakmatch_list[i]) :
                        peakmatch_list[i].append(current_newmatch)
                    else :
                        old_dSubScore = peakmatch_list[i][0][6]
                        if (dSubScore > old_dSubScore) :
                            peakmatch_list[i][0]=current_newmatch
                    bFindPeak = True
                    break
    return bFindPeak


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

def ClassifyBonds(current_mol, bBreakRing) :
    ring_bonds_list   = []
    linear_bonds_list = []
    iBondsNum = current_mol.GetNumBonds() 
    for i in range(iBondsNum) :
        current_bond = current_mol.GetBondWithIdx(i)
        if current_bond.IsInRing()  :
            if (bBreakRing) :
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
    #print len(current_fragments_list)
    if len(current_fragments_list) == 2 :
        bValidOperation = True
    else :
        bValidOperation = False
    #print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    return current_fragments_list, bValidOperation


def processKid(current_editable_mol, current_removebond_list, iCurrent_depth, allPeaks_list, peakmatch_list, sEnergy_Bond_dict, energy_denominator, bBreakRing) :
    #print current_editable_mol
    current_mol = current_editable_mol.GetMol()
    current_dMass, current_sFragmentFormula, current_smiles=DumpOneFragment(current_mol, current_removebond_list)
    bFindPeak = MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, current_removebond_list, sEnergy_Bond_dict, energy_denominator)
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(current_mol, bBreakRing)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []

    new_item = [[current_editable_mol, current_removebond_list, iCurrent_depth], current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid]
    return new_item

def CalculateEnergyDenominator(current_mol, sEnergy_Bond_dict) :
    energy_denominator = 0
    iBondsNum = current_mol.GetNumBonds()
    for i in range(iBondsNum) :
        current_bond = current_mol.GetBondWithIdx(i)
        sBeginAtom = current_bond.GetBeginAtom().GetSymbol()
        sEndAtom   = current_bond.GetEndAtom().GetSymbol()
        sType      = str(current_bond.GetBondType())
        if (sType == "SINGLE"):
            sRepresentType = "-"
        elif (sType == "DOUBLE"):
            sRepresentType = "="
        elif (sType == "TRIPLE"):
            sRepresentType = "~"
        #elif (sType == "AROMATIC"):
        #    sRepresentType = "*"
        else :
            sRepresentType = "("+sType+")"
        sCurrentBondType = sBeginAtom + sRepresentType + sEndAtom 
        if sCurrentBondType in sEnergy_Bond_dict :
            dCurrentBDE = sEnergy_Bond_dict[sCurrentBondType]
        else :
            dCurrentBDE = 348  #C-C
        energy_denominator += 1.0/dCurrentBDE
    if (energy_denominator == 0 ):
        energy_denominator = 1
    return energy_denominator

def ExhaustBonds(current_mol, allPeaks_list, peakmatch_list, sEnergy_Bond_dict, bBreakRing) :
    
    energy_denominator = CalculateEnergyDenominator(current_mol, sEnergy_Bond_dict)
    
    current_sFragmentFormula = AllChem.CalcMolFormula(current_mol)
    current_smiles = Chem.MolToSmiles(current_mol)
#    print Descriptors.ExactMolWt(current_mol), "NULL", AllChem.CalcMolFormula(current_mol)
    MapMass(Descriptors.ExactMolWt(current_mol), allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, [],sEnergy_Bond_dict, energy_denominator)
    iBondsNum = current_mol.GetNumBonds()
    #TreeLikeBreakBonds(current_mol, iBondsNum, allPeaks_list, 3, peakmatch_list)
    TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, 3, peakmatch_list, sEnergy_Bond_dict, energy_denominator, bBreakRing)


