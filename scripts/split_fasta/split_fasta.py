import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import reverse_complement, translate

inputFileName     = sys.argv[1]
outputFilePrefix  = sys.argv[2]
fileNum           = int(sys.argv[3])
seqNum            = int(sys.argv[4])

fileId   = 0
seqCount = 0
outputFile = open(outputFilePrefix+"_"+str(fileId)+".fasta", "w")

for seqRecord in SeqIO.parse( inputFileName   , "fasta" ) : 
	seqCount = seqCount + 1
	outputFile.write(">"+seqRecord.description+"\n")
	outputFile.write(seqRecord.seq.tostring()  +"\n")
	if ((seqCount == seqNum) and (fileId < (fileNum - 1))) :
		seqCount = 0
		fileId   = fileId + 1
		outputFile.close()
		outputFile = open(outputFilePrefix+"_"+str(fileId)+".fasta", "w")

outputFile.close()
