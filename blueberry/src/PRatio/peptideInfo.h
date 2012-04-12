
#ifndef PEPTIDEINFO_H
#define PEPTIDEINFO_H

#include <vector>
#include <iterator>
#include <algorithm>
#include "peptideRatio.h"
#include "peptideLabelFree.h"
#include "tinyxml.h"
#include "chromatogram.h"

using namespace std;

class PeptideRatio;
class PeptideLabelFree;

class PeptideInfo
{
	public:
		PeptideInfo();
		~PeptideInfo();
		
		// set all the values 
		void setValues( PeptideRatio * pPeptideRatio );
		void setValues( PeptideLabelFree * pPeptideLabelFree );

		// get and set functions for shared return variables between PeptideRatio and PeptideLabelFree
		void setFilename( string sFilenameInput );
		void setIdentifier( int iIdentifierInput );
		void setSequence( string sSequenceInput );
		void setChargeState( int iChargeStateInput );
		void setMaximumScore( float fScoreInput );
		void setValidity( bool bValidityInput );
		void setLocus( vector< string > vsLocusInput );

		string getFilename();
		int getIdentifier();
		string getSequence();
		int getChargeState();
		float getMaximumScore();
		bool getValidity();
		const vector< string > & getLocus();
		const vector< string > & getDescription();

		unsigned int getMS2Count() { return iMS2Count; };
		float getPeakTimeWidth() { return fPeakTimeWidth;};
		float getLeftValleyTime() { return fLeftValleyTime;};
		float getRightValleyTime() { return fRightValleyTime;};
		vector<float> getMS2Time() { return vfMS2Time;};
		vector<string> getAllIDfilename() { return vsAllIDfilename; };

		// get and set functions for return variables from PeptideRatio 
		void setPCALog2Ratio( double dPCALog2RatioInput );
		void setPCALog2SNR( double dPCALog2SNRInput );

		double getPCALog2Ratio();
		double getPCALog2SNR();

		// get and set functions for return variables from PeptideLabelFree
		double getPeakHeight(){return dPeakHeight;};
		double getPeakArea(){return dPeakArea;};
		double getPeakSNR(){return dPeakSNR;};
	private:

		// shared return variables between PeptideRatio and PeptideLabelFree
		string sFilename;
		int iIdentifier;
		string sSequence;
		int iChargeState;
		float fMaximumScore;
		bool bValidity;
		vector< string > vsLocus;
		vector< string > vsDescription;
		unsigned int iMS2Count;
		float fPeakTimeWidth;
		float fLeftValleyTime;
		float fRightValleyTime;
		vector<float> vfMS2Time;
		vector<string> vsAllIDfilename;

		// return variables from PeptideRatio 
		double dPCALog2Ratio;
		double dPCALog2SNR;

		// return variables from PeptideLabelFree
		double dPeakHeight;
		double dPeakArea;
		double dPeakSNR;

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
