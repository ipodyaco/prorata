#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math
import multiprocessing
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def parse_options (argv) :
    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                ["help", "input-database-filename", "output-database-filename"])
    sInput_Filename      = ""
    sOutput_Filename    = ""

    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-database-filename -o output-database-filename"
            sys.exit(0)
        if option in ("-i", "--input-database-filename"):
            sInput_Filename  = value
        if option in ("-o", "--output-database-filename"):
            sOutput_Filename = value
    
    if (sInput_Filename == "") :
        print "Please specify input file name"
        sys.exit(1)
    if (sOutput_Filename == "") :
        sOutput_Filename =  os.path.splitext(sInput_Filename)[0]+".neu.meb"
    
    return [sInput_Filename, sOutput_Filename]


def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]

    updated_smiles = Chem.MolToSmiles(mol,True)
    updated_mol    = Chem.MolFromSmiles(updated_smiles)

    if replaced:
        return (updated_mol, True)
    else:
        return (updated_mol, False)

def CalCompoundMass(sInchi) :
    current_mol  = Chem.MolFromInchi(sInchi)
    #reactions = InitiateNeutralisationReactions()
    neutral_mol, replaced  = NeutraliseCharges(Chem.MolToSmiles(current_mol))
    if (replaced) :
        current_mass = Descriptors.ExactMolWt(neutral_mol)
        sInchi       = Chem.MolToInchi(neutral_mol)
    else :
        current_mass = Descriptors.ExactMolWt(current_mol)
    return current_mass, sInchi, replaced

def ReadCompoundFile(compound_filename, output_filename) :
    compound_file = open(compound_filename)
    output_file   = open(output_filename, "w")
    bTitleLine = True
    bNewFormat = True
    for each_line in compound_file :
        each_line = each_line.strip()
        if ((each_line == "") or (each_line.startswith("#"))) :
            continue
        current_compound_info = each_line.split("\t")
        if ((len(current_compound_info) != 4) and (len(current_compound_info) != 5)) :
            print "illegal format", each_line
            sys.exit(1)
        elif (len(current_compound_info) == 5) :
            bNewFormat = False
        if (bTitleLine) and (bNewFormat):
            bTitleLine = False
            output_file.write(each_line+"\n")
            continue
            
        if (bNewFormat) :
            sInchi = current_compound_info[2]
        else :
            sInchi = current_compound_info[1]
        dMass, updated_inchi, bReplaced = CalCompoundMass(sInchi)
        #print dMass, updated_inchi
        if (bReplaced) :
            if (bNewFormat) :
                current_compound_info[2] = updated_inchi
            else :
                current_compound_info[1] = updated_inchi
                current_compound_info[2] = str(dMass)
            output_file.write(current_compound_info[0])
            for i in range(1, len(current_compound_info)) :
                output_file.write("\t"+current_compound_info[i])
            output_file.write("\n")
        else :
            output_file.write(each_line+"\n")

    compound_file.close()
    output_file.close()


def main(argv=None):
    if argv is None:
        argv = sys.argv
        [sInput_Filename, sOutput_Filename] = parse_options(argv)
        ReadCompoundFile(sInput_Filename, sOutput_Filename)
        
## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
