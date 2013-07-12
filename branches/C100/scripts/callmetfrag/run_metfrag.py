#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math


def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hw:r:o:t:a:m:",
                                    ["help",
                                     "working-dir",
                                     "realhit-filename",
				                     "output-dir",
                                     "thread-number",
                                     "ms-accuracy",
                                     "path-metfrag"])

    # Default working dir and config file
    working_dir   = ""
    output_dir    = ""
    realhit_filename   = ""
    path_metfrag  = ""
    iThreadNumber = 0
    dMsAccuracy   = -1.0
    path_metfrag  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w working-dir, -r realhit-filename  -o output-dir, -t thread-number, -a ms-accuracy, -m path-metfrag"
            sys.exit(0)
        if option in ("-w", "--working-dir"):
            working_dir  = value
            if (working_dir[-1] != os.sep) :
                working_dir += os.sep
        if option in ("-o", "--output-dir"):
            output_dir   = value
            if (output_dir[-1] != os.sep) :
                output_dir += os.sep
        if option in ("-r", "--realhit-filename") :
            realhit_filename  = value
        if option in ("-t", "--thread-number") :
            iThreadNumber = int(value)
        if option in ("-a", "--ms-accuracy") :
            dMsAccuracy = float(value)
        if option in ("-m", "--path-metfrag"):
            path_metfrag = value
    if ((working_dir=="")or(output_dir=="")or(realhit_filename=="")or(iThreadNumber<=0)or(dMsAccuracy<0.0)or(path_metfrag=="")) :
        print "please specify working-dir, output-dir, realhit-filename, ms-accuracy, and path-metfrag"
        sys.exit(1)

    return [working_dir, output_dir, realhit_filename, iThreadNumber, dMsAccuracy, path_metfrag]

def call_metfrag(working_dir, output_dir, iThreadNumber, dMsAccuracy, path_metfrag, corename, exact_mass) :
    command_str  = os.environ["JAVA_HOME"]+"/bin/java -jar "+path_metfrag+"  -d sdf  -L "+working_dir+corename+".sdf"
    command_str += " -D " + working_dir+corename+".mb  -S  "+corename+" -R "+output_dir+"  -n "+exact_mass
    command_str += " -a " + str(dMsAccuracy) + " -T " + str(iThreadNumber) + " -B "
    #print command_str
    os.system(command_str)

def RunTest (working_dir, output_dir, realhit_filename, iThreadNumber, dMsAccuracy, path_metfrag) :
    realhit_file = open(realhit_filename)
    for each_line in realhit_file :
        if (each_line.startswith("*\t")) :
            scan_info = each_line.split("\t")
            corename  = scan_info[1]
            exact_mass= scan_info[3]
            call_metfrag(working_dir, output_dir, iThreadNumber, dMsAccuracy, path_metfrag, corename, exact_mass)
    realhit_file.close()


def main(argv=None):
    if argv is None:
        argv = sys.argv
       		 # parse options
        [working_dir, output_dir, realhit_filename, iThreadNumber, dMsAccuracy, path_metfrag] = parse_options(argv)  	
        RunTest(working_dir, output_dir, realhit_filename, iThreadNumber, dMsAccuracy, path_metfrag) 

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()





