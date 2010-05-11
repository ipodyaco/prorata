#ifndef PROTEININDIRECTCOMPARISON_H
#define PROTEININDIRECTCOMPARISON_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "proteinDirectComparison.h"
#include "proteinInfo.h"
#include "proteinRatio.h"
#include "directoryStructure.h"
#include "proRataConfig.h"
#include "tinyxml.h"

using namespace std;

class ProteinIndirectComparison
{
	
	public:
		ProteinIndirectComparison();
		~ProteinIndirectComparison();
		
		bool runQuantification( string sInputLocus,
				ProteinDirectComparison directComparison0, ProteinDirectComparison directComparison1);
		void eraseProfileLikelihoodCurve();
		
		// set configurations
		void setName( string sNameConfig );
		void setNumerator( string sNumeratorConfig );
		void setDenominator( string sDenominatorConfig );	
		void setDirectComparisonName0( string sDirectComparisonName0Config);
		void setDirectComparisonName1( string sDirectComparisonName1Config );
		bool testDirectComparison( ProteinDirectComparison directComparison0, ProteinDirectComparison directComparison1);
		
		string getName();
		string getNumerator();
		string getDenominator();	
		string getDirectComparisonName0();
		string getDirectComparisonName1();	

		bool isValid();
		void getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef );
		
		// set functions for the configurations of indirect comparison
		// generally called by function readProRataCombineXML in ProjectInfo 
		static void setMaxCIwidth( double dMaxWidthCI)
		{ dMaxCIwidth = dMaxWidthCI; }		
		static void setLnLikehoodCutOffset( double dLnLikelihoodOffset )
		{ dLnLikehoodCutOffset = dLnLikelihoodOffset; }

	private:
		// static variables for the configurations of indirect comparison
		static double dMaxCIwidth;
		static double dLnLikehoodCutOffset;

		// the configurations for each indirect comparison
		string sName;
		string sNumerator;
		string sDenominator;
		// names of the two comparisons of this protein to be combined
		string sDirectComparisonName0;
		string sDirectComparisonName1;
		
		double dMLEMinLog2Ratio;
		double dMLEMaxLog2Ratio;
		double dLog2RatioDiscretization;
		
		/*
		 *  there are four types of combination
		 *  Indirect comparison = A:B
		 *  Common factor betwen the two direct comparison (DC) = X
		 *  1) indirect comparison =  DC0 - DC1     => { DC0 = A:X, DC1 = B:X } OR { DC0 = X:B, DC1 = X:A } 	
		 *  2) indirect comparison =  -(DC0 - DC1)  => { DC0 = X:A. DC1 = X:B } OR { DC0 = B:X, DC1 = A:X }	
		 *  3) indirect comparison =  DC0 + DC1     => { DC0 = A:X. DC1 = X:B } OR { DC0 = X:B, DC1 = A:X }	
		 *  4) indirect comparison =  -(DC0 + DC1)  => { DC0 = X:A. DC1 = B:X } OR { DC0 = B:X, DC1 = X:A }
		 *
		 *  iCombinationType = 0, if none of above is found in this indirect comparison	
		 */
		int iCombinationType;

		// return false, if iCombinationType = 0
		bool calCombinationType(string sDC0Numerator, string sDC0Denominator,
				string sDC1Numerator, string sDC1Denominator);

		double combineLog2Ratio( double dLog2RatioDC0, double dLog2RatioDC1 );


		// results for this protein's direct comparison		
		bool bValidity;
		string sLocus;
		double dLog2Ratio;
		double dLowerLimitCI;
		double dUpperLimitCI;

		// the profile likelihood curve		
		vector<double> vdLog2Ratio;
		vector<double> vdLnLikelihood;

};

#endif //PROTEININDIRECTCOMPARISON_H
