#ifndef FASTADATA_H
#define FASTADATA_H

#include <vector>
#include <cctype>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class FASTAdata
{

	public:
		FASTAdata();
		FASTAdata( string sFilename );
		~FASTAdata();

		// read a FASTA file and populate mLocus2Description and mLocus2Sequence
		bool readFASTA( string sFilename );

		// get the description for a locus
		string getDescription( string sLocus );

		// get the protein sequence for a locus
		string getProteinSequence( string sLocus );
		
		// compute the terminal distances for a peptide in a locus
		void computeDistances( string sLocus, string sPeptideSequence, int *piNterminalDistance, int *piCterminalDistance);

		void computeCoveragePlot( const string & sProteinSequence, string sPeptideSequence, int & iNtermCord, int & iCtermCord );
		
	private:
		// is this line a description line?
		bool isDesciptionLine( string sLine );

		// parse out locus and description and save them to mLocus2Description
		string processDescriptionLine(string sDescriptionLine);
		
		// if sSequenceLine is not a comment, then append it to sTempSequence
		string processSequenceLine( string sTempSequence, string sSequenceLine );
		
		// format the protein sequence and save it to mLocus2Sequence
		void saveSequence( string sLocus, string sSequence );

		// format the sequence: remove all non-alphabet letter and change all letter to upper case
		string formatSequence( string sSequence );
		
		map< string, string > mLocus2Description;
		map< string, string > mLocus2Sequence; 
		
};

#endif //FASTADATA_H
