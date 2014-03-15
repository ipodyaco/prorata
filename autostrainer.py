#!/usr/bin/python

import sys, getopt, warnings, os, re
import time
import itertools, copy, math
import pysam
from Bio import SeqIO

def parse_options(argv):

    opts, args = getopt.getopt(argv[1:], "hf:m:o:e:k:",
                                    ["help",
                             	     "input-fasta-file",
                                     "input-sam-file",
	                			     "output-file",
                                     "error-threshold",
                                     "variant-path-number"])

    output_file_path = ""
    iErrThreshold = 2
    iTopVariantPaths = 3

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-f input_fasta_file_name, -m input_sam_file_name -o output_file_path -e error threshold -k number of variant paths"
            sys.exit(1)
        if option in ("-f", "--input-fasta-file"):
            input_fasta = value
        if option in ("-m", "--input-same-file"):
            input_sam  = value
        if option in ("-o", "--output-file"):
            output_file_path = value
        if option in ("-e", "--error-threshold"):
            iErrThreshold = int(value)
        if option in ("-k", "--variant-path-number"):
            iTopVariantPaths = int (value)

    if ((input_sam == "" ) or (input_fasta == "")) :
        print "Please specify the fasta file and the sam file"
    if (output_file_path == ""):
        output_file_path = os.path.dirname(input_sam)
    if (not(output_file_path.endswith(os.sep))) :
        output_file_path += os.sep
    if (iErrThreshold <0 ):
        print "Error threshold can't be less than 0"
        sys.exit(1)
    if (iTopVariantPaths < 1) :
        print "Varint path number must be positive"
        sys.exit(1)


    return (input_fasta, input_sam, output_file_path, iErrThreshold, iTopVariantPaths)

def testsamfile(input_samfile) :
    seqnum = 0
    for alignedread in input_samfile.fetch():
        print alignedread
        print alignedread.pos, alignedread.qname,  alignedread.seq, alignedread.is_paired, alignedread.is_read1, alignedread.aligned_pairs
        print "**********************************"
        alignedread.seq ="yingfengwang"
        print alignedread
        seqnum +=1
        if (seqnum == 10) :
            break


def ReadContig(input_fasta_filename) :
    contig_record = SeqIO.parse(open(input_fasta_filename, "rU"), "fasta").next()
    contig_seq = str(contig_record.seq)
    contig_id  = str(contig_record.id)
    return contig_seq, contig_id

def pair_info(alignedread) :
    materead_name = ""
    read_name = alignedread.qname
    if (read_name == "") :
        print "no read name"
        sys.exit(1)

    if (alignedread.is_paired) :
            
        if alignedread.is_read1 :
            sReadId = "1"
            sMateReadId = "2"
        elif alignedread.is_read2 :
            sReadId = "2"
            sMateReadId = "1"
        else :
            print "incorrect pair info"
            sys.exit(1)

        if (read_name.endswith("/1")) :
            materead_name = read_name[:-1]+"2"
        elif (read_name.endswith("/2")) :
            materead_name = read_name[:-1]+"1"
        else :
            read_name += "/"+sReadId
            materead_name += "/"+sMateReadId
        if not(read_name.endswith(sReadId)) :
            print read_name
            print "incorrect pair info"
            sys.exit(1)
    #print read_name, materead_name
    return read_name, materead_name

def BuildRegionMap (input_samfile, sContigId, start_pos, end_pos) :
    region_list = [[] for i in range(end_pos-start_pos+1)]
    #print sContigId
    max_pos = 0
    for alignedread in input_samfile.fetch():
        #print alignedread
        #print alignedread.aligned_pairs
        align_info_list =  alignedread.aligned_pairs
        read_name = alignedread.qname
        read_seq = alignedread.seq
        for current_aligned_pair in align_info_list :
            iReadPos = current_aligned_pair[0]
            iRefPos  = current_aligned_pair[1]
            if (iRefPos > max_pos) :
                max_pos = iRefPos
            if ((iReadPos != None)  and (iRefPos != None)):
            # I have not observed case with iRefPos == None
                current_position_info = (read_seq[iReadPos:iReadPos+1], read_name, iReadPos)
            #    region_list[iRefPos-start_pos].append(current_position_info)
    print max_pos
    return region_list


def BuildReadDict (input_samfile) :
    read_dict = {}# all reads
    for alignedread in input_samfile.fetch():
        read_name, materead_name = pair_info(read_dict, alignedread)
        #read_seq_list = list (alignedread.seq)
        read_seq = alignedread.seq
        bIncluded = True
        current_read = [bIncluded, materead_name,  alignedread, read_seq]
        read_dict[read_name] = current_read

    return read_dict

def FindVariants(input_samfile, sRefSeq, iErrThreshold) :
    iRefLen = len(sRefSeq)
    ref_list = [[sRefSeq[i:i+1], 0, 0, 0, 0] for i in range(iRefLen)]
    for alignedread in input_samfile.fetch():
        align_info_list =  alignedread.aligned_pairs
        read_seq = alignedread.seq
        read_seq = read_seq.upper()
        for current_aligned_pair in align_info_list :
            iReadPos = current_aligned_pair[0]
            iRefPos  = current_aligned_pair[1]
            if (iRefPos >= iRefLen) :
                print iRefPos, "run out of ref index", iRefLen
                sys.exit(1)
            if ((iReadPos != None)  and (iRefPos != None)):
                current_NT = read_seq[iReadPos : iReadPos +1]
                if ((current_NT == "A") and (current_NT != ref_list[iRefPos][0])):
                    ref_list[iRefPos][1] += 1
                elif ((current_NT == "C") and (current_NT != ref_list[iRefPos][0])):
                    ref_list[iRefPos][2] += 1
                elif ((current_NT == "G") and (current_NT != ref_list[iRefPos][0])):
                    ref_list[iRefPos][3] += 1
                elif ((current_NT == "T") and (current_NT != ref_list[iRefPos][0])):
                    ref_list[iRefPos][4] += 1
    variants_list = []
    for i in range(iRefLen) :
        current_pos_info = ref_list[i]
        for j in range(1, 5) :
            if (current_pos_info[j] > iErrThreshold) :
                variants_list.append([i, current_pos_info[0], [0], [0], [0], [0], [0]])
                break
    print "Find", len(variants_list), "variant positions on total", iRefLen, "positions."
    del ref_list
    return variants_list

def ParseOutputFileName(input_sam_filename, output_path) :
    (input_sam_root, output_filename_ext) = os.path.splitext(input_sam_filename)
    input_sam_basename = os.path.basename(input_sam_root)
    output_filename_beg = output_path + input_sam_basename

    #print output_filename_beg, output_filename_ext
    return output_filename_beg, output_filename_ext

def OutputConsistentReads(variants_list, input_samfile, sRefSeq, iErrThreshold, output_filename_beg, output_filename_ext):
    ref_reads_filename = output_filename_beg + ".ref" + output_filename_ext
    ref_reads_file = pysam.Samfile(ref_reads_filename, "wh", template=input_samfile)
    NT_list = ["A", "C", "G", "T"]
    variant_read_dict = {}
    ref_variant_dict ={}
    iRefReadNum = 0
    for i in range(len(variants_list)):
        ref_variant_dict[variants_list[i][0]]= i
    for alignedread in input_samfile.fetch():
        bConsisitentRead = True
        align_info_list =  alignedread.aligned_pairs
        read_name = alignedread.qname
        read_seq = alignedread.seq
        for current_aligned_pair in align_info_list :
            iReadPos = current_aligned_pair[0]
            iRefPos  = current_aligned_pair[1]
            if ((iReadPos != None)  and (iRefPos != None)): 
                iVariantIndex = ref_variant_dict.get(iRefPos)
                if ( iVariantIndex != None) :
                    current_variant_info_list = variants_list[iVariantIndex]
                    sCurrentReadNt = read_seq[iReadPos:iReadPos+1]
                    updated_read_name, updated_materead_name = pair_info(alignedread)
                    if sCurrentReadNt not in NT_list : # not ACGT
                        sCurrentReadNt = variants_list[iVariantIndex][1]
                        read_seq = read_seq[:iReadPos] + sCurrentReadNt + read_seq[iReadPos+1:] 
                    if sCurrentReadNt in NT_list : 
                        iNtIndex = NT_list.index(sCurrentReadNt)
                    else : # reference NT not ACGT
                        iNtIndex = 4
                    current_variant_info_list[iNtIndex+2][0] += 1
                    current_variant_info_list[iNtIndex+2].append([updated_read_name, iReadPos, sCurrentReadNt])
                    if (bConsisitentRead) :
                        bConsisitentRead = False
                        bPicked4Path = False
                        variant_read_dict[updated_read_name] = [bPicked4Path, alignedread, read_seq, read_name, updated_materead_name]

        if (bConsisitentRead) :
            ref_reads_file.write(alignedread)
            iRefReadNum +=1

    ref_reads_file.close()
    print "The consistent sam file containing", iRefReadNum ,"reads has been done."
    #print len(variant_read_dict)
    return variant_read_dict


def PreprocessVariants(variants_list) :
# handle wild cards
    NT_list = ["A","C","G","T"]
    for i in range(len(variants_list)) :
    #    if (variants_list[i][0]== 21814) :
    #        print i, variants_list[i]
        current_pos_info = variants_list[i]
        iWildCardNumber  = current_pos_info[6][0]
        sRefGenomeNT     = current_pos_info[1]
        if (iWildCardNumber > 0):
            if sRefGenomeNT in NT_list : 
                # all wild cards go to the reference NT
                iNtIndex = NT_list.index(sRefGenomeNT)
                current_pos_info[2+iNtIndex][0] += iWildCardNumber
                # we don't update read seq in the read dict 
                for j in range(1, iWildCardNumber+1) :
                    current_pos_info[6][j][2] = sRefGenomeNT
                    current_pos_info[2+iNtIndex].append(current_pos_info[6][j])
                if (current_pos_info[2+iNtIndex][0] != (len(current_pos_info[2+iNtIndex]) -1)) :
                    print "wrong wild card case"
                    sys.exit(1)
            else:
                # all wild cards in the read go to the majority NT
                iNtNumber_list = [0,0,0,0]
                for j in range(4) :
                    iNtNumber_list[j] = current_pos_info[2+iNtIndex][0]
                iMajorityIndex = iNtNumber_list.index(max(iNtNumber_list))
                sMajorityNt    = NT_list[iMajorityIndex]
                current_pos_info[2+iMajorityIndex][0] += iWildCardNumber
                for j in range(1, iWildCardNumber+1) :
                    current_pos_info[6][j][2] = sMajorityNt
                    current_pos_info[2+iMajorityIndex].append(current_pos_info[6][j])
                if (current_pos_info[2+iMajorityIndex][0] != (len(current_pos_info[2+iMajorityIndex]) -1)) :
                    print "wrong wild card case"
                    sys.exit(1)

def RemovePickedReads (variants_list, picked_reads_dict, iNt_Choice_list) :
    for i in range(len(variants_list)) :
        current_read_list = variants_list[iNt_Choice_list[i]+2]
        for j in reversed(range(1, len(current_read_list))) :
            updated_read_name = current_read_list[j][0]
            if ( updated_read_name  in  picked_reads_dict ) :
                del current_read_list[j]


def ClassifyReads(current_variant_position_list, rejected_reads_dict, variant_read_dict) :
    Reads_list = [[], [], [], []]
    for i in range(4) :
        ntbased_read_list = current_variant_position_list[i+2][1:]
        for j in range(len(ntbased_read_list)) :
            current_read_name = ntbased_read_list[j][0]
            if ((current_read_name not in rejected_reads_dict) and (current_read_name in variant_read_dict)):
                read_info_list = variant_read_dict[current_read_name]
                current_mate_read_name = read_info_list[4]
                if ((current_mate_read_name == "") or (current_mate_read_name not in variant_read_dict)) :
                    read_info_list[4] = ""
                    Reads_list[i].append(current_read_name)
                elif (current_mate_read_name not in rejected_reads_dict) :
                    Reads_list[i].append(current_read_name)
#                else :
#                    rejected_reads_dict[current_read_name] = True

    return Reads_list


def PickVariants (Reads_list, rejected_reads_dict, picked_reads_dict, variant_read_dict) :
    iVariantNum_list = [len(Reads_list[i]) for i in range(4)]
    iPicked_id = iVariantNum_list.index(max(iVariantNum_list))
    if (iVariantNum_list[iPicked_id] == 0) :
        return False
    else:# update rejected_reads_dict and  picked_reads_dict
        for i in range(4) :
            for j in range(iVariantNum_list[i]) :
                current_read_name = Reads_list[i][j]
                read_info_list = variant_read_dict[current_read_name]
                current_mate_read_name = read_info_list[4]
                if (i == iPicked_id):
                    picked_reads_dict[current_read_name]   = True
                else :
                    rejected_reads_dict[current_read_name] = True
                    if current_read_name in picked_reads_dict :
                        del picked_reads_dict[current_read_name]
                if (current_mate_read_name != "") :
                    if current_mate_read_name not in variant_read_dict :
                        read_info_list[4] = ""
                    else :
                        if (i == iPicked_id):
                            picked_reads_dict[current_mate_read_name]   = True
                        else :
                            rejected_reads_dict[current_read_name] = True
                            if current_mate_read_name in picked_reads_dict :
                                del picked_reads_dict[current_mate_read_name]


    return True

def BuildVariantPaths(variants_list, variant_read_dict, variant_path_id, iNt_Choice_list, picked_reads_dict, rejected_reads_dict) :
    if (variant_path_id == 0) :
        PreprocessVariants(variants_list)
    else :
        RemovePickedReads (variants_list, picked_reads_dict, iNt_Choice_list)
    for i in range(len(variants_list)):
        current_variant_position_list = variants_list[i]
        Reads_list = ClassifyReads(current_variant_position_list, rejected_reads_dict, variant_read_dict)
        #if (i==352):
        #    print Reads_list
        bEnoughRead = PickVariants (Reads_list, rejected_reads_dict, picked_reads_dict, variant_read_dict)
    #    if not(bEnoughRead) :
    #        print "Not enough reads for generating variant path", variant_path_id+1, "on position", variants_list[i][0], "(index starting from 0) with totally", len(variant_read_dict), "reads left."
            #print i, len(variant_read_dict)
    #        sys.exit(1)
    # update variants_read_dict
    for updated_read_name in picked_reads_dict.keys() :
        current_read_info_list = variant_read_dict[updated_read_name]
        current_read_info_list[0] = True

    print len(picked_reads_dict), "reads for variant path", variant_path_id + 1 
    


def OutputVariantPaths(variant_read_dict, variant_path_id, input_samfile, output_filename_beg, output_filename_ext) :
    variant_path_filename = output_filename_beg + ".var" + str(variant_path_id+1)  + output_filename_ext
    variant_path_file = pysam.Samfile(variant_path_filename, "wh", template=input_samfile)
    for current_updated_readname, current_readinfo_list in variant_read_dict.items() :
        bPicked4Path = current_readinfo_list[0]
        if (bPicked4Path) :
            alignedread = current_readinfo_list[1]
            variant_path_file.write(alignedread)
            del  variant_read_dict[current_updated_readname]

    variant_path_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
       		 # parse options
        input_fasta_filename, input_sam_filename, output_path, iErrThreshold, iTopVariantPaths = parse_options(argv)  
        #iErrThreshold = 2
        #iTopVariantPaths = 3
        output_filename_beg, output_filename_ext = ParseOutputFileName(input_sam_filename, output_path) 
        input_samfile = pysam.Samfile( input_sam_filename, "r" )
        sRefSeq, sContigId = ReadContig(input_fasta_filename)
        sRefSeq = sRefSeq.upper()
        variants_list = FindVariants(input_samfile, sRefSeq, iErrThreshold)
        variant_read_dict = OutputConsistentReads(variants_list, input_samfile, sRefSeq, iErrThreshold, output_filename_beg, output_filename_ext)
        iNt_Choice_list = [-10 for i in range(len(variants_list))]
        

        for i in range(iTopVariantPaths) :
            picked_reads_dict   = {}
            rejected_reads_dict = {}
            BuildVariantPaths(variants_list, variant_read_dict, i, iNt_Choice_list, picked_reads_dict, rejected_reads_dict)
            OutputVariantPaths(variant_read_dict, i, input_samfile, output_filename_beg, output_filename_ext)
        
        #BuildReadDict(input_samfile)
        #BuildRegionMap(input_samfile, sContigId, 0, 5)

        #print sRefSeq 
        #testsamfile(input_samfile)
        #test2col(input_sam, ref_genomes)
        input_samfile.close()

## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
