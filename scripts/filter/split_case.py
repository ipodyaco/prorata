#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math, random


input_filename = sys.argv[1]
output_path    = sys.argv[2]
each_subfile_case_num = int(sys.argv[3])

current_case_num = 0
file_id          = 0

input_file = open(input_filename)

case_buff = ""

for each_line in input_file :
    each_line = each_line.strip()
    case_buff += each_line+"\n"
    if (each_line == "//") :
        current_case_num += 1
    if (current_case_num == each_subfile_case_num) :
        output_file = open(output_path+str(file_id)+".txt", "w")
        output_file.write(case_buff)
        output_file.close()
        current_case_num = 0
        file_id += 1
        case_buff = ""

output_file = open(output_path+str(file_id)+".txt", "w")
output_file.write(case_buff)
output_file.close()

input_file.close()


