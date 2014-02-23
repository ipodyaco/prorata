#!/usr/bin/python
import sys, getopt, warnings, os, re

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                                     "input-file",
                                     "output-dir"])

    sInputFileName   = ""
    sOutputDir       = "" 

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input_file_name -o output_directory"
            sys.exit(1)
        if option in ("-i", "--input-file"):
            sInputFileName = value
        if option in ("-o", "--output-dir"):
            sOutputDir = value
            if sOutputDir[-1] != os.sep :
                sOutputDir += os.sep

    if (sInputFileName == "") or (sOutputDir == "") :
        print "please specify input file name and output directory"
        sys.exit(1)

    return sInputFileName, sOutputDir


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
        die("Program exit!")

    return file_list


def getColumnId(sColumnNameLine, ColumnName) :
# This function returns column id (starting with 0) according to column name
# sColumnNameLine: the comlumn title line; ColumnName: column name
	ColumnName_list = sColumnNameLine.split("\t")
	try:
		iColumnId = ColumnName_list.index(ColumnName)
	except ValueError:
		print "can't find column "+ColumnName
		sys.exit(0)
	#print iColumnId
	return iColumnId

def CallRDKit(sSmiles, sTempFileName) :
    
    current_mol = Chem.MolFromSmiles(sSmiles)
        #current_mol = Chem.MolFromSmiles(sSmiles, sanitize=False)
    try:
        AllChem.Compute2DCoords(current_mol)
    except:
        print sSmiles
        current_mol = Chem.MolFromSmiles(sSmiles, sanitize=False)
        AllChem.Compute2DCoords(current_mol)

    Draw.MolToFile(current_mol, sTempFileName, size=(50, 50))


def DrawFragments(sCurrentSmiles_list, sTempFragmentFigureName_list, sOutputDir) :
    iFragmentId = 0
    for each_smiles in sCurrentSmiles_list :
        #print each_smiles
        sTempFileName = sOutputDir+"fragment_"+str(iFragmentId)+".png"
        sTempFragmentFigureName_list.append(sTempFileName)
        CallRDKit(each_smiles, sTempFileName)
        iFragmentId += 1
       # current_mol = Chem.MolFromSmiles(each_smiles)
       # AllChem.Compute2DCoords(current_mol)
       # Draw.MolToFile(current_mol, sTempFileName, size=(50, 50))

def DrawOneDiagram(sCurrentIdentifier, dCurrentMofZ_Identified, dCurrentIntensity_Identified, dCurrentMofZ_unIndentified, dCurrentIntensity_unIndentified, sCurrentSmiles_list, sOutputDir) :
    sTempFragmentFigureName_list = []
    DrawFragments(sCurrentSmiles_list, sTempFragmentFigureName_list, sOutputDir)

    x_max = max([max(dCurrentMofZ_Identified) , max(dCurrentMofZ_unIndentified)])
    y_max = max([max(dCurrentIntensity_Identified),max(dCurrentIntensity_unIndentified)])

    x_min = min([min(dCurrentMofZ_Identified) , min(dCurrentMofZ_unIndentified)])
    
#    print dCurrentMofZ_Identified
#    print dCurrentIntensity_Identified

    plt.bar(dCurrentMofZ_Identified, dCurrentIntensity_Identified, width=0.2, color='red', edgecolor='none')
    plt.bar(dCurrentMofZ_unIndentified, dCurrentIntensity_unIndentified, width=0.2, color='black', edgecolor='none')

    plt.xlabel('m/z')
    plt.ylabel('relative intensity')

    plt.axis([x_min-50, x_max+10, 0, y_max+0.3])
    plt.draw()
    #plt.show()
    plt.savefig(sOutputDir+sCurrentIdentifier+".svg", format="svg")
    for each_tempfile in sTempFragmentFigureName_list :
        os.remove(each_tempfile)

def DrawAllDiagrams(sInputFileName, sOutputDir) :
    input_file = open(sInputFileName)
    iLineNum = 0
    for each_line in input_file:
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        iLineNum += 1
        if (iLineNum == 1) :
            iIdentifier_ColumnId = getColumnId(each_line, "Identifier") + 1
        elif (iLineNum == 2) :
            iMofZ_ColumnId = getColumnId(each_line, "m/z")
            iIntensity_ColumnId = getColumnId(each_line, "Intensity")
            iSmiles_ColumnId = getColumnId(each_line, "SMILES")
        elif (each_line.startswith("M\t")) :
            if (iLineNum > 3) :
                DrawOneDiagram(sCurrentIdentifier, dCurrentMofZ_Identified, dCurrentIntensity_Identified, dCurrentMofZ_unIndentified, dCurrentIntensity_unIndentified, sCurrentSmiles_list, sOutputDir)
            currentMebInfo_list = each_line.split("\t")
            sCurrentIdentifier = currentMebInfo_list[iIdentifier_ColumnId]
            dCurrentMofZ_Identified = []
            dCurrentIntensity_Identified = []
            dCurrentMofZ_unIndentified = []
            dCurrentIntensity_unIndentified = []
            sCurrentSmiles_list = []
        else :
            currentPeakInfo_list = each_line.split("\t")
            sCurrentCompound = currentPeakInfo_list[iSmiles_ColumnId]
            dCurrentMofZ = float(currentPeakInfo_list[iMofZ_ColumnId])
            dCurrentIntensity = float(currentPeakInfo_list[iIntensity_ColumnId])
            if (sCurrentCompound == "NA") :
                dCurrentMofZ_unIndentified.append(dCurrentMofZ)
                dCurrentIntensity_unIndentified.append(dCurrentIntensity)
            else :
                dCurrentMofZ_Identified.append(dCurrentMofZ)
                dCurrentIntensity_Identified.append(dCurrentIntensity)
                sCurrentSmiles_list.append(sCurrentCompound)
    if (iLineNum > 3) :
        DrawOneDiagram(sCurrentIdentifier, dCurrentMofZ_Identified, dCurrentIntensity_Identified, dCurrentMofZ_unIndentified, dCurrentIntensity_unIndentified, sCurrentSmiles_list, sOutputDir)

    input_file.close()
    


def main(argv=None):
    if argv is None:
        argv = sys.argv
        sInputFileName, sOutputDir = parse_options(argv)
        DrawAllDiagrams(sInputFileName, sOutputDir)

if __name__ == "__main__":
    main()



