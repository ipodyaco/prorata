#!/usr/bin/python

import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator
from Bio import SeqIO
from Bio.Seq import Seq

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
    matrix_filename_list = get_file_list_with_ext(working_dir, ".matrix.txt")
    return [matrix_filename_list, output_dir]



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

def parseEmatrixLine(Ematrix_line) :
    line_info_list = Ematrix_line.split("\t")
    current_Ematrix_row = []
    for i in range(2, len(line_info_list)) :
        current_field = line_info_list[i]
        current_equation = current_field.split("(")[0]
        equation_info_list = current_equation.split("=")
        if (int(equation_info_list[1]) == 1) :
            current_Ematrix_row.append(int(equation_info_list[0]))
    return current_Ematrix_row

def parseSmatrixLine(Smatrix_line, i_column_number) :
    current_Smatrix_row = [0] * i_column_number
    line_info_list = Smatrix_line.split("\t")
    for i in range(2, len(line_info_list)) :
        current_field = line_info_list[i]
        current_equation = current_field.split("(")[0]
        equation_info_list = current_equation.split("=")
        current_Smatrix_row[int(equation_info_list[0])] = float(equation_info_list[1])
    return current_Smatrix_row


def ParseMatrixFile(matrix_filename) :
    isEmatrix = True
    Ematrix_list = []
    Smatrix_list = []
    i_column_number = 0
    i_row_number    = 0
    i_matrix_number = 0
    matrix_file = open(matrix_filename)
    for each_line in matrix_file:
        each_line = each_line.strip()
        if ((each_line.startswith("#")) or (each_line == "") ) :
            continue
        if each_line.startswith("+") :
            i_matrix_number += 1
            if (i_matrix_number > 2) :
                break
            if (i_column_number != i_row_number) :
                print "column number ", i_column_number, " doesn't match row nubmer", i_row_number
                sys.exit(1)
            i_column_number = 0
            i_row_number    = 0
            if (isEmatrix) :
                title_info_list = each_line.split("\t")
                s_matrix_name= title_info_list[1]
                s_matrix_description  = title_info_list[2]
                pos = s_matrix_name.find("_")
                s_model_name = s_matrix_name[pos+1:]
                isEmatrix = False
        if each_line.startswith("@") :
            i_column_number += 1
        if each_line.startswith("*") :
            i_row_number += 1
            if (i_matrix_number == 1) :
                current_Ematrix_row = parseEmatrixLine(each_line)
                Ematrix_list.append(current_Ematrix_row)
            if (i_matrix_number == 2) :
                current_Smatrix_row = parseSmatrixLine(each_line, i_column_number)
                Smatrix_list.append(current_Smatrix_row)
    matrix_file.close()
#    print Ematrix_list
#    print Smatrix_list
    return [Ematrix_list, Smatrix_list, s_matrix_name, s_matrix_description]

def PathSearch(Ematrix_list, Smatrix_list, current_path_list, end_vertex_id, d_current_path_weight) :
    next_vertex_list = Ematrix_list[current_path_list[-1]]
    if (next_vertex_list == []) :
        if (current_path_list[-1] == end_vertex_id) :
            PathSet_list.append(current_path_list)
            PathWeight_list.append(d_current_path_weight)
            PathLength_list.append(len(current_path_list))
  #          Path_Vertex_Covered_list.append([False]*len(current_path_list))
    else:
        for each_next_vertex in next_vertex_list :
            d_weight_addition = 0
            for each_path_vertex in current_path_list :
               d_weight_addition += Smatrix_list[each_next_vertex][each_path_vertex]
               d_weight_addition += Smatrix_list[each_path_vertex][each_next_vertex]
            PathSearch(Ematrix_list, Smatrix_list, current_path_list+[each_next_vertex], end_vertex_id, d_current_path_weight+d_weight_addition)


def PathEnumeration(Ematrix_list, Smatrix_list) :
    global PathSet_list
    global PathWeight_list
    global PathLength_list
 #   global Path_Vertex_Covered_list # vertices covered by other selected pathed, True, covered
    PathSet_list = []
    PathWeight_list = []
    PathLength_list = []
 #   Path_Vertex_Covered_list = []
    PathSearch(Ematrix_list, Smatrix_list, [0], len(Ematrix_list)-1, 0)

def NextSet(Cover_Set_list) :
    best_path_id = -100
    best_average_weight = -100
    for i in range(len(PathSet_list)):
        if i not in Cover_Set_list :
            if (PathLength_list[i] > 0) :
                current_average_weight = PathWeight_list[i]/PathLength_list[i]
                if (current_average_weight > best_average_weight) :
                    best_path_id = i
                    best_average_weight = current_average_weight
    if (best_path_id < 0) :
        print "fail to pick a path"
        sys.exit(0)
    return best_path_id

def GreedySetCover(i_vertex_number) :
    Cover_Set_list = []
    Cover_Vertex_list = [False] * i_vertex_number
    while False in Cover_Vertex_list :
        if (len(Cover_Set_list) == len(PathSet_list) ):
            print "can't cover all vertices"
            sys.exit(1)
        Selected_Set_id  = NextSet(Cover_Set_list)
        Cover_Set_list.append(Selected_Set_id)
        for i in range(len(PathSet_list[Selected_Set_id])) :
            current_vertex_id = PathSet_list[Selected_Set_id][i]
            if (Cover_Vertex_list[current_vertex_id] == False ) :
                Cover_Vertex_list[current_vertex_id] = True
                for j in range(len(PathSet_list)) :
                    if j not in Cover_Set_list :
                        if current_vertex_id in PathSet_list[j] :
                            PathLength_list[j] -= 1

        
    return Cover_Set_list


def SearchSetCover(Ematrix_list, Smatrix_list, matrix_filename, output_dir, s_matrix_name, s_matrix_description) :
    (matrix_filename_head, matrix_filename_tail) = os.path.split(matrix_filename)
    (matrix_filename_root, matrix_filename_ext)  = os.path.splitext(matrix_filename_tail)
    (matrix_filename_root, matrix_filename_ext)  = os.path.splitext(matrix_filename_root)
    output_filename = output_dir + matrix_filename_root + ".pathset.txt"
    output_file = open(output_filename, "w")
    
    PathEnumeration(Ematrix_list, Smatrix_list)
    #PathSet_ori_list = PathSet_list
    Cover_Set_list =  GreedySetCover(len(Ematrix_list))

    output_file.write("+\tPathSet_"+s_matrix_name+"\t"+s_matrix_description+"\n")
    for i in range(len(Ematrix_list)) :
        output_file.write("@\t"+str(i)+"\tvertex_"+str(i)+"\n")
    for i in range(len(Cover_Set_list)) :
        output_file.write("*\tpath_"+str(i))
        for each_vertex_id in PathSet_list[Cover_Set_list[i]] :
            output_file.write("\t"+str(each_vertex_id))
        output_file.write("\n")

    output_file.close()
#    print PathSet_list, PathWeight_list, PathLength_list, Path_Vertex_Covered_list


## +------+
## | Main |
## +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
        [matrix_filename_list, output_dir] = parse_options(argv)

    for current_matrix_filename in matrix_filename_list :
        [Ematrix_list, Smatrix_list, s_matrix_name, s_matrix_description] = ParseMatrixFile(current_matrix_filename)
        SearchSetCover(Ematrix_list, Smatrix_list, current_matrix_filename, output_dir, s_matrix_name, s_matrix_description)



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()


