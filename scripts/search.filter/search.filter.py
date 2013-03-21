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




## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [parse_filename_list, output_dir] = parse_options(argv)

    for each_search_filename in search_filename_list :
        handleSearchFile(each_search_filename, output_dir) 

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()




