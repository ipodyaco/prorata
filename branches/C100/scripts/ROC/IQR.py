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
        sOutputFileName = os.path.splitext(sInputFileName)[0] + ".iqr"
    
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

def IQRcalculation(dCurrentScore_list) :
    dCurrentScore_list.sort()
    iItemNumber = len(dCurrentScore_list)
#    if (iItemNumber > 20) :
#        dCurrentScore_list = dCurrentScore_list[-20:]
#        iItemNumber = 20
    if (iItemNumber < 4) :
        dIQR_1st = -100
        dIQR_2nd = -100
        dIQR_3rd = -100
    else :
        iFirstIndex = int(iItemNumber * 0.25) 
        iThirdIndex = int(iItemNumber * 0.75)
        dFirst = dCurrentScore_list[iFirstIndex]
        dThird = dCurrentScore_list[iThirdIndex]
        #print iItemNumber, iFirstIndex, iThirdIndex
        dHighest = dCurrentScore_list[-1]
        d2ndHighest = dCurrentScore_list[-2]
        d3rdHighest = dCurrentScore_list[-3]
        #print dHighest, d2ndHighest
        dDiffOne = dHighest - dThird
        dDiffTwo = d2ndHighest -dThird
        dDiffThree = d3rdHighest -dThird
        dIQRDistance = dThird - dFirst
        if (dIQRDistance == 0) :
            if dDiffOne > 0 :
                dIQR_1st = 100000
            else :
                dIQR_1st =-100
            if dDiffTwo > 0 :
                dIQR_2nd = 100000
            else :
                dIQR_2nd = -100
            if dDiffThree > 0 :
                dIQR_3rd = 100000
            else :
                dIQR_3rd = -100
        else :
            dIQR_1st = dDiffOne/dIQRDistance
            dIQR_2nd = dDiffTwo/dIQRDistance
            dIQR_3rd = dDiffThree/dIQRDistance
    return dIQR_1st, dIQR_2nd, dIQR_3rd

def ReadScore(sInputFileName) :
    input_file = open(sInputFileName)
    iCaseNum = 0
    bRealHit = False
    dRealHit_IQRs_list = []
    dWrongHit_IQRs_list= []
    dCurrentScore_list = []
    for each_line in input_file :
        each_line = each_line.strip()
        if (each_line == ""):
            continue
        if (each_line.startswith(">")):
            if (dCurrentScore_list) :
                dIQR_1st, dIQR_2nd, dIQR_3rd = IQRcalculation(dCurrentScore_list)
                if (bRealHit) :
                    #print dIQR_1st, dIQR_2nd, dIQR_3rd, dIQR_1st-dIQR_2nd, dIQR_2nd-dIQR_3rd
                    if (len(dCurrentScore_list) > 2) :
                        dCurrentScore_list.sort()
                        dRealHit_IQRs_list.append(dIQR_1st)
                        dWrongHit_IQRs_list.append(dIQR_2nd)
               # print dIQR, bRealHit
            dCurrentScore_list = []
            iCaseNum += 1
            bRankOne = True
            if each_line.startswith(">1\t"):
                bRealHit = True
            else :
                bRealHit = False
            if (len(each_line.split("\t")) > 2) :
                score_idx = 2
            else :
                score_idx = 1
        else :
            dScore = float(each_line.split("\t")[score_idx])
            dCurrentScore_list.append(dScore)
    if (bRealHit) :
        if (len(dCurrentScore_list) > 2) :
            dCurrentScore_list.sort()
            dRealHit_IQRs_list.append(dIQR_1st)
            dWrongHit_IQRs_list.append(dIQR_2nd)
    input_file.close()
    print iCaseNum
    print len(dRealHit_IQRs_list)
    return dRealHit_IQRs_list, dWrongHit_IQRs_list, iCaseNum


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
        dRealHit_IQRs_list, dWrongHit_IQRs_list, iCaseNum = ReadScore(sInputFileName)
        curve_point = CalculateFDR(dRealHit_IQRs_list, dWrongHit_IQRs_list)
        WriteOutput(curve_point, sOutputFileName)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()



