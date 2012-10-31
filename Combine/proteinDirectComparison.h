#ifndef PROTEINDIRECTCOMPARISON_H
#define PROTEINDIRECTCOMPARISON_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <iomanip>
#include <fstream>

#include "proteinReplicate.h"
#include "proteinInfo.h"
#include "proteinRatio.h"
#include "proRataConfig.h"
#include "tinyxml.h"

using namespace std;


class ProteinDirectComparison
{
	
	public:
		ProteinDirectComparison();
		~ProteinDirectComparison();

		bool runQuantification( string sInputLocus );
		void eraseProfileLikelihoodCurve();
		bool getProfileLikelihoodCurve(vector<double> & vdLog2RatioRef, vector<double> & vdLnLikelihoodRef );
		const vector< ProteinReplicate > & getProteinReplicateVector();

		int getQuantifiedPeptides();

		void setName( string sNameConfig );
		void setNumerator( string sNumeratorConfig );
		void setDenominator( string sDenominatorConfig );

		string getName();
		string getNumerator();
		string getDenominator();

	//	double getMLEMinLog2Ratio();
	//	double getMLEMaxLog2Ratio();
		double getLog2RatioDiscretization();
		
		string getLocus();
		bool isValid();
		void getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef );
		int getValidReplicateNumber();
		int getValidPeptideNumber();
		
		
		bool addReplicate( ProteinReplicate tempProteinReplicate );

		// set functions for the configurations of direct comparison
		// generally called by function readProRataCombineXML in ProjectInfo 
		static void setValidReplicateCICutoff( int iValidCINumberCutoff )
		{iValidReplicateCICutoff = iValidCINumberCutoff; }
		static void setMinimumPeptideNumber( int iNumberPeptide )
		{ iMinimumPeptideNumber  = iNumberPeptide; }
		static void setMaxCIwidth( double dMaxWidthCI)
		{ dMaxCIwidth = dMaxWidthCI; }
		static void setLnLikehoodCutOffset( double dLnLikelihoodOffset )
		{ dLnLikehoodCutOffset = dLnLikelihoodOffset; }
		static void setRequireCIoverlapped( bool bRequireCIoverlappedInput )
		{ bRequireCIoverlapped = bRequireCIoverlappedInput; }
		static void setMaxLogRatioDifference( double MaxLogRatioDifferenceInput )
		{ dMaxLogRatioDifference = MaxLogRatioDifferenceInput; }

	private:
		// static variables for the configurations of direct comparison
		static int iValidReplicateCICutoff;
		static int iMinimumPeptideNumber;
		static double dMaxCIwidth;
		static double dLnLikehoodCutOffset;
		static bool bRequireCIoverlapped;
		static double dMaxLogRatioDifference;

		// the configurations read from ProRata_Combine.XML
		string sName;
		string sNumerator;
		string sDenominator;

		// CONFIG from replicate QPR
		double dMLEMinLog2Ratio;
		double dMLEMaxLog2Ratio;
		double dLog2RatioDiscretization;
		
		// the replicates of this protein to be combined
		vector< ProteinReplicate > vReplicate;

		// results for this protein's direct comparison
		string sLocus;
		bool bValidity;
		double dLog2Ratio;
		double dLowerLimitCI;
		double dUpperLimitCI;
		int iQuantifiedPeptides;
		int iValidReplicateNumber;
		int iValidPeptideNumber;

		// the profile likelihood curve
		vector<double> vdLog2Ratio;
		vector<double> vdLnLikelihood;

		bool bIsCIOverlapped(double dCI0lower, double dCI0upper, double dCI1lower, double dCI1upper);
};

#endif //PROTEINDIRECTCOMPARISON_H
