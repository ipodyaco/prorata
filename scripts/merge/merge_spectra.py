#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
import time
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
            print "-i input-file, -o output-file "
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-file"):
            output_filename = value

    if (input_filename == "") or (output_filename == ""):
        print "Please specify -i or -o"
        sys.exit(1)

    return (input_filename, output_filename)

def ReadInputFile(input_filename) :
    input_file = open(input_filename)
    all_spectra_list =[]
    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith("*")):
            current_spectrum = []
            current_peaks_list = []
            current_inchi = each_line.split("\t")[4].strip()
            current_spectrum.append(current_inchi)
            current_spectrum.append(each_line)
        if (each_line.startswith("+")) :
            current_mofz      = float(each_line[1:].split(" ")[0].strip())
            current_intensity = float(each_line[1:].split(" ")[-1].strip())
            current_peaks_list.append([current_mofz, current_intensity])
        if (each_line.startswith("//")) :
            current_spectrum.append(current_peaks_list)
            all_spectra_list.append(current_spectrum)
    input_file.close()
    return all_spectra_list

def MergeTwoSpectra(spectrum_A, spectrum_B, merge_range) :
    spectrum_1st = spectrum_A[:]
    spectrum_2nd = spectrum_B[:]
    spectrum_1st.sort(key=lambda m:m[0])
    spectrum_2nd.sort(key=lambda m:m[0])
    for each_peak in spectrum_1st :
        for i in range(len(spectrum_2nd))[::-1] :
            if (math.fabs(each_peak[0] - spectrum_2nd[i][0]) < merge_range) :
                each_peak[0] = (each_peak[0] + spectrum_2nd[i][0])/2.0
                each_peak[1] = max(each_peak[1], spectrum_2nd[i][1])
                del spectrum_2nd[i]
    combined_list = spectrum_1st + spectrum_2nd 
    combined_list.sort(key=lambda m:m[0])
    return combined_list

def MergeSpectra(all_spectra_list, merge_range) :
    all_combined_spectra_list = []
    all_inchi_list = []
    for each_spectrum in all_spectra_list :
        current_inchi = each_spectrum[0]
        if current_inchi not in all_inchi_list :
            all_inchi_list.append(current_inchi)
            all_combined_spectra_list.append(each_spectrum)
        else :
            current_id = all_inchi_list.index(current_inchi)
            all_combined_spectra_list[current_id][2] = MergeTwoSpectra(all_combined_spectra_list[current_id][2],each_spectrum[2], merge_range)
    return all_combined_spectra_list

def WriteOutput(all_combined_spectra_list, output_filename) :
    output_file = open(output_filename, "w")
    for each_spectrum in all_combined_spectra_list :
        output_file.write(each_spectrum[1]+"\n")
        for each_peak in each_spectrum[2] :
            output_file.write("+"+str(each_peak[0])+" "+str(each_peak[1])+"\n")
        output_file.write("//\n")

    output_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        input_filename ,  output_filename = parse_options(argv)  
        merge_range = 0.01
        all_spectra_list = ReadInputFile(input_filename)
        #print all_spectra_list[0], all_spectra_list[-1], len(all_spectra_list)
        all_combined_spectra_list = MergeSpectra(all_spectra_list, merge_range)
        #print all_combined_spectra_list[0], all_combined_spectra_list[-1], len(all_combined_spectra_list
        WriteOutput(all_combined_spectra_list, output_filename)

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()

