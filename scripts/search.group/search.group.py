#!/usr/bin/python

import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator


## Parse options
def parse_options(argv):
    opts, args = getopt.getopt(argv[1:], "hw:o:f:n:", ["help", "working-dir", "output-dir", "read-filename", "sub-sequence-number"])
    # Default working dir and config file
    working_dir = "./"
    output_dir  = ""
    read_filename = ""
    sub_sequence_number = 0
    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-o outputdirectory -w workingdirectory -f read_filename -n subset sequence_number"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-dir"):
            output_dir = value 
        if option in ("-f", "--read-filename"):
            read_filename = value
        if option in ("-n", "--sub-sequence-number"):
            sub_sequence_number = int(value)
    # only -w is provided
    if (output_dir == "") :
        output_dir = working_dir       
    if (sub_sequence_number == 0) or (read_filename == "") :
        print "please specify read_filename and subset sequence_number"
        sys.exit(1)
    filter_filename_list = get_file_list_with_ext(working_dir, ".search.filter.txt")
    return [filter_filename_list, read_filename, output_dir, sub_sequence_number]



## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):
    file_list = []
    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):
            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)
        file_list = sorted(file_list)
    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        sys.exit(1)  
    return file_list

def identifyAllColumnID(first_line) :
    column_names = first_line.split("\t")
    for i in range(len(column_names)) :
        if (column_names[i] == "Pfam_ID") :
            PfamID_columnID = i
        elif (column_names[i] == "Peptide_ID") :
            PeptideID_columnID = i
        elif (column_names[i] == "score") :
            score_columnID = i
        elif (column_names[i] == "hmmfrom") :
            hmmfrom_columnID = i
        elif (column_names[i] == "hmm_to") :
            hmmto_columnID = i
        elif (column_names[i] == "alifrom") :
            alifrom_columnID = i
        elif (column_names[i] == "ali_to") :
            alito_columnID = i
        elif (column_names[i] == "Model_alignment") :
            Modelalignment_columnID = i
        elif (column_names[i] == "Peptide_alignment") :
            Peptidealignment_columnID = i
    return [PfamID_columnID, PeptideID_columnID, score_columnID, hmmfrom_columnID, hmmto_columnID, alifrom_columnID, alito_columnID, Modelalignment_columnID, Peptidealignment_columnID] 


## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [filter_filename_list, read_filename, output_dir, sub_sequence_number] = parse_options(argv)

    for each_parse_filename in parse_filename_list :
        handleParseFile(each_parse_filename, output_dir) 

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
