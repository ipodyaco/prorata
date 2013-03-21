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
    search_filename_list = get_file_list_with_ext(working_dir, ".search.output.txt")
    return [search_filename_list, output_dir]



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

def IdentifyScore (infoLine) :
    featureCount = 0
    infoLine.strip(" ")
    allFeature = infoLine.split(" ")
    for i in range(len(allFeature)) :
        if (len(allFeature[i]) > 0 ) :
            if (allFeature[i][:1].isdigit()) :
                featureCount = featureCount + 1
                if (featureCount == 2) :
                    Score = allFeature[i]
                elif (featureCount == 4) :
                    cEvalue = allFeature[i]
                elif (featureCount == 5) :
                    iEvalue = allFeature[i]
                elif (featureCount == 6) :
                    hmmFrom = allFeature[i]    
                elif (featureCount == 7) :
                    hmmTo   = allFeature[i]
                elif (featureCount == 8) :
                    aliFrom = allFeature[i]
                elif (featureCount == 9) :
                    aliTo   = allFeature[i]
                elif (featureCount == 10) :
                    envFrom = allFeature[i]
                elif (featureCount == 11) :
                    envTo   = allFeature[i]
                    break
    return [Score, cEvalue, iEvalue, hmmFrom, hmmTo, aliFrom, aliTo, envFrom, envTo]


def handleSearchFile(search_filename, output_dir) :
    (search_filename_head, search_filename_tail)= os.path.split(search_filename)
    (search_filename_root, search_filename_ext) = os.path.splitext(search_filename_tail)
    (search_filename_root, search_filename_ext) = os.path.splitext(search_filename_root)
    output_filename = output_dir + search_filename_root + ".parse.txt"
    #print output_filename
    search_file = open(search_filename)
    output_file = open(output_filename, "w")

    output_file.write("Pfam_ID\tPeptide_ID\tscore\tc-Evalue\ti-Evalue\thmmfrom\thmm_to\talifrom\tali_to\tenvfrom\tenv_to\tModel_alignment\tPeptide_alignment\n")
    
    isDomain = False
    lineCount= -1

    for eachline in search_file :
        eachline = eachline.strip()
        if (eachline == "") :
            continue
        if (eachline.startswith("Query:")) :
            isDomain = False
            queryInfo = eachline.split(" ")
            modelId  = queryInfo[7]
        if (eachline.startswith(">>")) : # a hit
            peptideId = eachline[2:].strip()
            isDomain = True
            lineCount= 0
        if (isDomain) :
            lineCount = lineCount + 1
        if (lineCount == 4) : # score, range, and other info
            scoreInfo = IdentifyScore (eachline)
        if (isDomain) : 
            if (eachline.startswith(modelId)):
                aliInfo     = eachline.split(" ")
                modelAlignment = aliInfo[-2]
            if (eachline.startswith(peptideId)):
                aliInfo     = eachline.split(" ")
                peptideAlignment = aliInfo[-2]
                isDomain = False
                output_file.write(modelId+"\t"+peptideId+"\t"+scoreInfo[0]+"\t"+scoreInfo[1]+"\t"+scoreInfo[2]+"\t"+scoreInfo[3]+"\t")
                output_file.write(scoreInfo[4]+"\t"+scoreInfo[5]+"\t"+scoreInfo[6]+"\t"+scoreInfo[7]+"\t"+scoreInfo[8]+"\t")
                output_file.write(modelAlignment+"\t"+peptideAlignment+"\n")

    search_file.close()
    output_file.close()


## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [search_filename_list, output_dir] = parse_options(argv)

    for each_search_filename in search_filename_list :
        handleSearchFile(each_search_filename, output_dir) 

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



