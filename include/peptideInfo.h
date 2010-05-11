
#ifndef PEPTIDEINFO_H
#define PEPTIDEINFO_H

#include <vector>
#include <iterator>
#include <algorithm>
#include "peptideRatio.h"
#include "tinyxml.h"
#include "chromatogram.h"

using namespace std;

class PeptideRatio;

class PeptideInfo
{
	public:
		PeptideInfo();
		~PeptideInfo();
		
		bool setPeptideRatio( PeptideRatio * pPeptideRatio );

		void setValues( PeptideRatio * pPeptideRatio );
		void setFilename( string sFilenameInput );
		void setIdentifier( int iIdentifierInput );
		void setSequence( string sSequenceInput );
		void setChargeState( int iChargeStateInput );
		void setMaximumScore( float fScoreInput );
		void setValidity( bool bValidityInput );
		void setPCALog2Ratio( double dPCALog2RatioInput );
		void setPCALog2SNR( double dPCALog2SNRInput );
		void setLocus( vector< string > vsLocusInput );

		string getFilename();
		int getIdentifier();
		string getSequence();
		int getChargeState();
		float getMaximumScore();
		bool getValidity();
		double getPCALog2Ratio();
		double getPCALog2SNR();
		const vector< string > & getLocus();
		const vector< string > & getDescription();
		
	private:
		string sFilename;
		int iIdentifier;
		string sSequence;
		int iChargeState;
		float fMaximumScore;
		bool bValidity;
		double dPCALog2Ratio;
		double dPCALog2SNR;
		vector< string > vsLocus;
		vector< string > vsDescription;
};

class LessPeptideInfo
{
	public:
		LessPeptideInfo();
		LessPeptideInfo( string sKeyInput ) { sKey = sKeyInput; }
		
		void setKey( string sKeyInput ) { sKey = sKeyInput; }
		string getKey() { return sKey; }
		
		/*
		 * the PeptideInfo pointers can be sorted by a number of member variables:
		 * "sequence"			pPeptide1->getSequence()
		 * "log2Ratio"			pPeptide1->getPCALog2Ratio()
		 * "log2SN"			pPeptide1->getPCALog2SNR()
		 */
		bool operator() ( PeptideInfo * pPeptide1, PeptideInfo * pPeptide2 ) const;
		
	private:
		string sKey;
};


#endif //PEPTIDEINFO_H
