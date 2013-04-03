#!/usr/bin/python

import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator
from Bio import SeqIO
from Bio.Seq import Seq


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
            print "-o outputdirectory -w workingdirectory -f read_filename -n subset_sequence_number"
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


def sort_filenames(filename_list) :
    all_filename_list = []
    for each_filename in filename_list :
        (each_filename_head, each_filename_tail) = os.path.split(each_filename)
        each_filename_root = each_filename_head.split(".")[0]
        each_filename_index= int(each_filename_root.split("_")[-1])
        all_filename_list.append([each_filename_index, each_filename])
    all_filename_list.sort(key = lambda x: x[0])
    return [each_filename_combine[1] for each_filename_combine in all_filename_list]


def AddRecordToFile(feature_list, current_PfamID, curent_read_seq, output_dir) :
    output_filename = output_dir+current_PfamID+".model.group.txt"
    if (os.path.exists(output_filename)) :
        output_file = open(output_filename, "w")
        output_file.write("Pfam_ID\tPeptide_ID\tscore\thmmfrom\thmm_to\talifrom\tali_to\tModel_alignment\tPeptide_alignment\tRead_Seq\n")
    else :
        output_file = open(output_filename, "a")
    output_file.write(feature_list[0]+"\t"+feature_list[1]+"\t"+feature_list[2]+"\t"+feature_list[3]+"\t"+feature_list[4]+"\t")
    output_file.write(feature_list[5]+"\t"+feature_list[6]+"\t"+feature_list[7]+"\t"+feature_list[8]+"\t"+currentReadSeq+"\n")
    output_file.close()

def handleFilterFile(current_read_list, current_filter_filename, output_dir) :
    
    is_title_line = True
    filter_file = open(current_filter_filename)

    for each_line in filter_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if (is_title_line) :
            [PfamID_columnID, PeptideID_columnID, score_columnID, hmmfrom_columnID, hmmto_columnID, alifrom_columnID, alito_columnID, 
                Modelalignment_columnID, Peptidealignment_columnID] = identifyAllColumnID(filter_file)
            is_title_line = False
        else :
            feature_list = each_line.split("\t")
            current_PfamID = feature_list[PfamID_columnID]
            current_PeptideID = feature_list[PeptideID_columnID]
            current_readID    = current_PeptideID[:-3]
            curent_read_seq   = ""
            for each_read in current_read_list :
                if (each_read.description == current_readID ) :
                   current_read_seq = each_read.seq.tostring()
                   break
            if (current_read_seq == ""):
                print "can't find seq: "+ current_readID
                sys.exit(1)
            else :
                AddRecordToFile(feature_list, current_PfamID, curent_read_seq, output_dir)

    filter_file.close()



## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [filter_filename_list, read_filename, output_dir, sub_sequence_number] = parse_options(argv)

    filter_filename_list = sort_filenames(filter_filename_list)
    read_iterator = SeqIO.parse( read_filename, "fasta" ) 
    current_read_list = []
    i_read_count = 0
    i_filter_filename_count = 0
    for each_read in read_iterator :
        current_read_list.append(each_read)
        i_read_count += 1
        if (i_read_count == sub_sequence_number) :
            current_filter_filename  = filter_filename_list[i_filter_filename_count]
            i_filter_filename_count += 1
            handleFilterFile(current_read_list, current_filter_filename, output_dir)
            i_read_count = 0
            current_read_list = []
    if (len(current_read_list) > 0) :
        handleFilterFile(current_read_list, filter_filename_list[i_filter_filename_count], output_dir)

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
