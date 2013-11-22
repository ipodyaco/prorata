#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <string>

void testf()
{
    seqan::FaiIndex faiIndex;
    //string sContigFileName = "test/case2/GWA2_AR10_Chromosome_circular_complete-FINAL_recirc.fasta";
    std::string sContigFileName = "test/case1/FINAL_GWA2_AR10_MAPPED_READS.fasta";
    std::string sSamFileName    = "test/case1/1.sam";
/*    
    seqan::BamStream bamStreamIn(sSamFileName.c_str());
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    bamStreamOut.header = bamStreamIn.header;
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        readRecord(record, bamStreamIn);
        writeRecord(bamStreamOut, record);
    }
    */
    
    if (read(faiIndex, sContigFileName.c_str()) != 0)
	if (build(faiIndex, sContigFileName.c_str()) != 0)
	    std::cerr<<"error"<<std::endl; //Error
    seqan::BamStream bamIO(sSamFileName.c_str());
    seqan::String<unsigned> mapping;
    resize(mapping, length(bamIO.header.sequenceInfos), 0);
    
    for (unsigned i = 0; i < length(bamIO.header.sequenceInfos); ++i)
    {
        seqan::CharString seqName = bamIO.header.sequenceInfos[i].i1;
	std::cout<<bamIO.header.sequenceInfos[i].i2<<std::endl;
        if (!getIdByName(faiIndex, seqName, mapping[i]))
        {
            std::cerr << "ERROR: Sequene "
                      << bamIO.header.sequenceInfos[i].i1
                      << "unknown in FASTA Index.\n";
            //return 1; 
        }
    }    
    
}

int main(int argc, char **argv) 
{
    
    std::cout << "Hello, world!" << std::endl;
    testf();
    return 0;
}
