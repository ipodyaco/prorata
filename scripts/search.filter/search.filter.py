#!/usr/bin/python

import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator


## Parse options
def parse_options(argv):
    opts, args = getopt.getopt(argv[1:], "hw:o:", ["help", "working-dir", "output-dir"])
    # Default working dir and config file
    working_dir = "./"
    output_dir  = ""
    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-o outputdirectory -w workingdirectory"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-dir"):
            output_dir = value   
    # only -w is provided
    if (output_dir == "") :
        output_dir = working_dir                
    parse_filename_list = get_file_list_with_ext(working_dir, ".search.parse.txt")
    return [parse_filename_list, output_dir]



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
        elif (column_names[i] == "c-Evalue") :
            cEvalue_columnID = i
        elif (column_names[i] == "i-Evalue") :
            iEvalue_columnID = i
        elif (column_names[i] == "hmmfrom") :
            hmmfrom_columnID = i
        elif (column_names[i] == "hmm_to") :
            hmmto_columnID = i
        elif (column_names[i] == "alifrom") :
            alifrom_columnID = i
        elif (column_names[i] == "ali_to") :
            alito_columnID = i
        elif (column_names[i] == "envfrom") :
            envfrom_columnID = i
        elif (column_names[i] == "env_to") :
            envto_columnID = i
        elif (column_names[i] == "Model_alignment") :
            Modelalignment_columnID = i
        elif (column_names[i] == "Peptide_alignment") :
            Peptidealignment_columnID = i
    return [PfamID_columnID, PeptideID_columnID, score_columnID, cEvalue_columnID, iEvalue_columnID, hmmfrom_columnID, hmmto_columnID, alifrom_columnID, alito_columnID, envfrom_columnID, envto_columnID, Modelalignment_columnID, Peptidealignment_columnID] 


def handleParseFile(parse_filename, output_dir) :
    (parse_filename_head, parse_filename_tail)= os.path.split(parse_filename)
    (parse_filename_root, parse_filename_ext) = os.path.splitext(parse_filename_tail)
    (parse_filename_root, parse_filename_ext) = os.path.splitext(parse_filename_root)
    output_filename = output_dir + parse_filename_root + ".filter.txt"
    input_file  = open(parse_filename)
    output_file = open(output_filename, "w")

    is_first_line = True

    readID_list = []
    score_list = []
    best_hit_list=[]

    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (is_first_line) :
            is_first_line = False
            [PfamID_columnID, PeptideID_columnID, score_columnID, cEvalue_columnID, iEvalue_columnID, hmmfrom_columnID, hmmto_columnID, alifrom_columnID, alito_columnID, envfrom_columnID, envto_columnID, Modelalignment_columnID, Peptidealignment_columnID] = identifyAllColumnID(each_line) 
            output_file.write("Pfam_ID\tPeptide_ID\tscore\thmmfrom\thmm_to\talifrom\tali_to\tModel_alignment\tPeptide_alignment\n")
        else :
            hit_info = each_line.split("\t")
            current_hit = hit_info[PfamID_columnID]+"\t"+hit_info[PeptideID_columnID]+"\t"+hit_info[score_columnID]+"\t"
            current_hit = current_hit+hit_info[hmmfrom_columnID]+"\t"+hit_info[hmmto_columnID]+"\t"+hit_info[alifrom_columnID]+"\t"
            current_hit = current_hit+hit_info[alito_columnID]+"\t"+hit_info[Modelalignment_columnID]+"\t"
            current_hit = current_hit+hit_info[Peptidealignment_columnID]
            current_readID = hit_info[PeptideID_columnID][:-3]
            current_score  = float(hit_info[score_columnID])
            if (current_readID in readID_list):
                current_index = readID_list.index(current_readID)
                if (current_score > score_list[current_index]) :
                    score_list[current_index]    = current_score
                    best_hit_list[current_index] = current_hit
            else :
                readID_list.append(current_readID)
                score_list.append(current_score)
                best_hit_list.append(current_hit)              
    for each_hit in best_hit_list :
        output_file.write(each_hit+"\n")

    input_file.close()
    output_file.close()


## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [parse_filename_list, output_dir] = parse_options(argv)

    for each_parse_filename in parse_filename_list :
        handleParseFile(each_parse_filename, output_dir) 

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()




