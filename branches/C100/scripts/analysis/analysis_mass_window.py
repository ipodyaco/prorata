#!/usr/bin/python

import sys, getopt, warnings, os, re
import time
import itertools, copy, math

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

def ReadOtherInfo(input_filename) :

    input_file =open(input_filename)

    all_info_list = []

    bEligibleEntry = False
    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (each_line.startswith(">1\t")):
            bEligibleEntry = True
        elif (bEligibleEntry) :
            bEligibleEntry = False
            current_info = each_line.split("\t")[-1]
            all_info_list.append(current_info)
    input_file.close()

    return all_info_list

def OutputMassWindowInfo(output_filename, all_info_list, mass_window) :
    
    output_file = open(output_filename, "w")
    for each_info in all_info_list :
        current_distribution = [0 for x in range(len(mass_window))]
        sOutputStr = ""
        sDistributionInfo = each_info.split(";")[1]
        if (sDistributionInfo != "") :  
            sDistribution_Window_list = sDistributionInfo.split(",")
            for each_window in sDistribution_Window_list :
                window_info_list = each_window.split(":")
                window_index = mass_window.index(int(window_info_list[0]))
                current_distribution[window_index] = int(window_info_list[1])
        for each_count in current_distribution :
            sOutputStr += str(each_count)+"\t"
        output_file.write(sOutputStr[:-1]+"\n")

    output_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        input_filename ,  output_filename = parse_options(argv)  
        all_info_list = ReadOtherInfo(input_filename)
        OutputMassWindowInfo(output_filename, all_info_list,[-1,0,1,2])

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
