#!/usr/bin/python

import sys, getopt, warnings, os, re
import itertools, copy, math



def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                                     "input-file",
				                     "output-file"])


    # Default working dir and config file
    sInputFileName   = ""
    sOutputFileName  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-file -o output-file"
            sys.exit(0)
        if option in ("-i", "--input-file"):
            sInputFileName  = value
        if option in ("-o", "--output-dir"):
            sOutputFileName = value

    if ( sInputFileName  == ""  )  :
        print "please specify input filename"
        sys.exit(1)

    if ( sOutputFileName == "") :
        sOutputFileName = os.path.splitext(sInputFileName)[0] + ".fdr"
    
    return [sInputFileName, sOutputFileName]

## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        sys.exit(1)

    return file_list

def parseInchi (swholeinchi) :
    inchi_info = swholeinchi.split("/")
    #if (len(inchi_info) < 3) :
    #    revalue = ""
    #else :
    revalue = inchi_info[1]
    return swholeinchi

def getColumnId(sColumnNameLine, ColumnName) :
    ColumnName_list = sColumnNameLine.split("\t")
    try:
        iColumnId = ColumnName_list.index(ColumnName)
    except ValueError:
        print "can't find column "+ColumnName
        sys.exit(0)
    return iColumnId


def ReadScore(sInputFileName):
    input_file = open(sInputFileName)
    iCaseNum = 0
    bTitleLine = True
    dRealHit_Scores_list = []
    dWrongHit_Scores_list= []

    sCurrentScanInfo = ""
    sFirstInchi_sub      = ""
    bSecond = True

    for each_line in input_file :
        each_line = each_line.strip()
        if ((each_line == "") or (each_line.startswith("#"))) :
            continue
        if (bTitleLine) :
            bTitleLine = False
            iFilenameColumnId = getColumnId(each_line, "Filename")
            iScanNumberColumnId = getColumnId(each_line, "ScanNumber")
            iRankColumnId  = getColumnId(each_line, "Rank")
            iScoreColumnId = getColumnId(each_line, "Score")
            iIdentifiedInChIColumnId = getColumnId(each_line, "IdentifiedInChI")
            iExplainedPeaksColumId = getColumnId(each_line, "ExplainedPeaks")
            continue
        hit_info_list = each_line.split("\t")
        sCurrentFilename = hit_info_list[iFilenameColumnId]
        sCurrentScanNumber = hit_info_list[iScanNumberColumnId]
        iCurrentRank = int(hit_info_list[iRankColumnId])
       # print iCurrentRank, sCurrentScanNumber
        dCurrentScore = float(hit_info_list[iScoreColumnId])
        sCurrentInchi = hit_info_list[iIdentifiedInChIColumnId]
        sCurrentExplainedPeaks = hit_info_list[iExplainedPeaksColumId]

        sLineScanInfo = sCurrentScanNumber+"@"+sCurrentFilename
        sLineInchi_sub= parseInchi(sCurrentInchi)
        iLineIdentifiedPeaks = int(sCurrentExplainedPeaks.split("/")[0])
        if (sLineScanInfo == sCurrentScanInfo) :
            if (sFirstInchi_sub != sLineInchi_sub) :
                print sFirstInchi_sub
            if ((sFirstInchi_sub != sLineInchi_sub) and (iCurrentRank > 1)  and (bSecond)):
                dWrongHit_Scores_list.append(dCurrentScore/iFirstIdentifiedPeaks)
                bSecond = False
                iCaseNum += 1
        else :
            sCurrentScanInfo = sLineScanInfo
            iFirstIdentifiedPeaks = iLineIdentifiedPeaks
            if (iFirstIdentifiedPeaks == 0) :
                iFirstIdentifiedPeaks = 1
            if (iCurrentRank != 1) :
                print "wrong rank", sCurrentScanInfo
            bSecond = True
            iCaseNum += 1
            sFirstInchi_sub = sLineInchi_sub
            dRealHit_Scores_list.append(dCurrentScore/iFirstIdentifiedPeaks)

    input_file.close()
    return dRealHit_Scores_list, dWrongHit_Scores_list, iCaseNum

def ReadScore_old(sInputFileName) :
    input_file = open(sInputFileName)
    iCaseNum = 0
    bRankOne = False
    bRankTwo = False
    bRealHit = False
    dRealHit_Scores_list = []
    dWrongHit_Scores_list= []
    
    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == ""):
            continue
        if (each_line.startswith(">")):
            iCaseNum += 1
            bRankOne = True
            if each_line.startswith(">1\t"):
                bRealHit = True
            else :
                bRealHit = False
            if (len(each_line.split("\t")) > 2) :
                inchi_idx = 3
                score_idx = 2
            else :
                inchi_idx = 4
                score_idx = 1
        elif (bRankOne) :
            dCurrentScore = float(each_line.split("\t")[score_idx])
           # print each_line
            sFirstInchi   = each_line.split("\t")[inchi_idx]
            sFirstInchi_sub = parseInchi(sFirstInchi)
            bRankOne = False
            bRankTwo = True
            if (bRealHit) :
#                bRealHit = False
                dRealHit_Scores_list.append(dCurrentScore)
#                print dCurrentScore, each_line.split("\t")[-1]
            #else :
            #    dWrongHit_Scores_list.append(dCurrentScore)
   #         print dCurrentScore
        elif (bRankTwo) :
            bRankTwo = False
            if (bRealHit) :
                dCurrentScore = float(each_line.split("\t")[score_idx])
                sSecondInchi  = each_line.split("\t")[inchi_idx]
                sSecondInchi_sub = parseInchi(sSecondInchi)
                if (sSecondInchi_sub != sFirstInchi_sub) :
                    dWrongHit_Scores_list.append(dCurrentScore)
                else :
                    bRankTwo = True
    input_file.close()
   # print iCaseNum
    return dRealHit_Scores_list, dWrongHit_Scores_list, iCaseNum


def WriteOutput(curve_point, sOutputFileName) :
    output_file = open(sOutputFileName, "w")
    for each_point in curve_point :
        output_file.write(str(each_point[0])+" "+str(each_point[1])+" "+str(each_point[2])+" "+str(each_point[3])+"\n")

    output_file.close()

def CalculateFDR(dRealHit_Scores_list, dWrongHit_Scores_list):
    d_allscores_list = dRealHit_Scores_list + dWrongHit_Scores_list
    d_allscores_list.append(-100)
    d_allscores_list.sort()
    curve_point = []
    for each_threshold in d_allscores_list :
        iGoodHit = 0
        iBadHit  = 0
        for each_score in dRealHit_Scores_list :
            if each_score >= each_threshold :
                iGoodHit += 1
        for each_score in dWrongHit_Scores_list :
            if each_score >= each_threshold :
                iBadHit  += 1
        currentFDR = iBadHit/float(iGoodHit+iBadHit)
        curve_point.append([each_threshold, currentFDR, iGoodHit, iBadHit])
    return curve_point

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        [sInputFileName, sOutputFileName] = parse_options(argv)      
        dRealHit_Scores_list, dWrongHit_Scores_list, iCaseNum = ReadScore(sInputFileName)
        curve_point = CalculateFDR(dRealHit_Scores_list, dWrongHit_Scores_list)
        WriteOutput(curve_point, sOutputFileName)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



