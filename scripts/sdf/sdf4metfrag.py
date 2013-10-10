#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hd:r:o:",
                                    ["help",
                                     "db-filename",
                                     "realhit-filename",
				                     "output-dir"])

    # Default working dir and config file
    database_filename  = ""
    output_dir    = ""
    realhit_filename   = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-d database_filename, -r realhit-filename  -o output-dir"
            sys.exit(0)
        if option in ("-d", "--db-filename"):
            database_filename  = value
        if option in ("-o", "--output-dir"):
            output_dir   = value
        if option in ("-r", "--realhit-filename"):
            realhit_filename  = value

    if ((database_filename == "") or (output_dir == "") or (realhit_filename == "")) :
        print "please specify database file name, output dir, and realhit filename"
        sys.exit(1)

    return [database_filename, output_dir, realhit_filename]

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


def Generate_mb_file(output_dir, realhit_filename) :
    mass_list = []
    bIsMBFileOpen = False
    realhit_file = open(realhit_filename)
    
    for each_line in realhit_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith("*")) :
            hit_info = each_line.split("\t")
            hit_name = hit_info[1]
            hit_mass = hit_info[3]
            hit_inchi= hit_info[4]
            mass_list.append([hit_name, hit_mass, hit_inchi])
            if (bIsMBFileOpen) :
                current_mb_file.write(peak_str[:-1])
                current_mb_file.close()
            current_mb_file = open(output_dir+os.sep+hit_name+".mb", "w")
            bIsMBFileOpen = True
            peak_str = ""
        if (each_line.startswith("+")) :
            peak_str += each_line[1:]+"\n"
            
    if (bIsMBFileOpen) :
        current_mb_file.write(peak_str[:-1])
        current_mb_file.close()
    realhit_file.close()
    return mass_list

def WriteToSdf(current_sdf_file, compound_internal_id, compound_inchi) :
    current_mol =  Chem.MolFromInchi(compound_inchi)
    current_mol.SetProp("UniqueID", compound_internal_id)
    current_mol.SetProp("Inchi", compound_inchi)
    current_sdf_file.write(current_mol)

def Generate_sdf(database_filename, output_dir, mass_list) :
    precursor_accuracy = 0.01
    for each_entry in mass_list :
        current_name = each_entry[0]
        current_mass = float(each_entry[1]) # the proton mass have been removed by another script
        current_inchi= each_entry[2]
        current_sdf_file = Chem.SDWriter(output_dir + os.sep + current_name + ".sdf")
        Compound_list    = ReadCompoundFile(database_filename)
        QueryCompound_list = GetRelatedCompound(Compound_list, current_mass, precursor_accuracy)
        for each_compound in QueryCompound_list :
            compound_internal_id = each_compound[0]
            compound_inchi = each_compound[1]
            WriteToSdf(current_sdf_file, compound_internal_id, compound_inchi)
        current_sdf_file.close()

def main(argv=None):
    if argv is None:
        argv = sys.argv
       		 # parse options
        [database_filename, output_dir, realhit_filename] = parse_options(argv)  	
        mass_list = Generate_mb_file(output_dir, realhit_filename) 
        Generate_sdf(database_filename, output_dir, mass_list)

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
