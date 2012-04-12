
#ifndef PROTEININFO_H
#define PROTEININFO_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>

#include "peptideRatio.h"
#include "peptideInfo.h"
#include "proteinRatio.h"

using namespace std;

class ProteinRatio;

class ProteinInfo
{
	public:
		ProteinInfo();
		~ProteinInfo();
		
		void computeLabelFree();

		bool setProteinRatio( ProteinRatio * pProteinRatioInput );

		void addPeptideInfo( PeptideInfo * pPeptideInfoInput );
		
		void setLocus( string sLocusInput );
		void setDescription( string sDescriptionInput);
		void setLog2Ratio( double dLog2RatioInput);
		void setLowerLimitCI( double dLowerLimitCIInput);
		void setUpperLimitCI( double dUpperLimitCIInput );
		
		void setValidPeptides( int  iValidPeptidesInput );
		void setQuantifiedPeptides( int iQuantifiedPeptidesInput );
		void setIdentifiedPeptides( int iIdentifiedPeptidesInput );

		vector< PeptideInfo * > & getPeptideInfo();
		vector< string > getPeptideSequence();
		string getLocus();
		string getDescription();
		double getLog2Ratio();
		double getLowerLimitCI();
		double getUpperLimitCI();
		int getValidPeptides();

		// this overloaded function is used in proteinReplicate
		int getValidPeptides(double dLowerLimitCI, double dUpperLimitCI);
		int getQuantifiedPeptides();
		int getIdentifiedPeptides();
		bool getValidity();

		// accessors for label-free variables
		double getTotalPeakHeight(){ return dTotalPeakHeight;};
		double getTotalPeakArea(){ return dTotalPeakArea;};
		int getMS2SpectralCounts(){ return iMS2SpectralCounts;};

	private:
		vector< PeptideInfo * > vpPeptideInfo;
		bool bValidity;
		string sLocus;
		string sDescription;

		// label free variables
		double dTotalPeakHeight;
		double dTotalPeakArea;
		int iMS2SpectralCounts;

		// metablic labeling variables
		double dLog2Ratio;
		double dLowerLimitCI;
		double dUpperLimitCI;
		
		// number of identified peptides 
		int iIdentifiedPeptides;
		
		// number of quantified peptides
		// a quantified peptide is an identified peptide with a log-S/N above the cutoff
		int iQuantifiedPeptides;
		
		// number of valid peptides
		// a valid peptide is a quantified peptide with a log-ratio within protein confidence interval
		int iValidPeptides;
};

class LessProteinInfo
{
	public:
		LessProteinInfo();
		LessProteinInfo( string sKeyInput ) { sKey = sKeyInput; }
		
		void setKey( string sKeyInput ) { sKey = sKeyInput; }
		string getKey() { return sKey; }
		
		/*
		 * the proteinInfo pointers can be sorted by a number of member variables:
		 * the sKey specifies which variables to be used to sort
		 * the accepted values for sKay and the corresponding accessors of the variables are
		 * "locus"			pProtein1->getLocus()
		 * "log2Ratio"			pProtein1->getLog2Ratio()
		 * "widthCI"			pProtein1->getUpperLimitCI() - pProtein1->getLowerLimitCI()
		 * "description"		pProtein1->getDescription()
		 * "peptidesQuantified"		pProtein1->getQuantifiedPeptides()
		 */
		bool operator() ( ProteinInfo * pProtein1, ProteinInfo * pProtein2 ) const;
		
	private:
		string sKey;
};

#endif //PROTEININFO_H
