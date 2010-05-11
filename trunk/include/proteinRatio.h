
#ifndef PROTEINRATIO_H
#define PROTEINRATIO_H

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "proRataConfig.h"

using namespace std;

// constants
const double PNORM_COEFF1 = 0.5;
const double PNORM_COEFF2 = 2.0;
const double BASE = 2.0;
// const double LOGLIKEOFFSET = 1.92;

class ProteinRatio
{
	public:
		ProteinRatio();
		~ProteinRatio();

		// get the inputs and compute the results
		bool process( const vector<double> & vdPeptideLog2RatioInput, 
				const vector<double> & vdPeptideLog2SNInput );

		// accessors for the results
		double getProteinLog2Ratio()
		{ return dProteinLog2Ratio; }

		double getLowerLimitCI()
		{ return dLowerLimitCI; }

		double getUpperLimitCI()
		{ return dUpperLimitCI; }

		const vector< double > & getPeptideLog2Ratio()
		{ return vdPeptideLog2Ratio; }
		
		const vector< double > & getPeptideLog2SN()
		{ return vdPeptideLog2SN; }

		const vector< double > & getLog2RatioBin()
		{ return vdLog2RatioBin; }
			
		const vector< double > & getLnLikelihood()
		{ return vdLnLikelihood; }
		
		double getMaxLnLikelihood()
		{ return dMaxLnLikelihood; }
		
		double getLnLikelihoodCutoff()
		{ return dLnLikelihoodCutoff; }

		int getLnLikelihoodSize()
		{ return vdLnLikelihood.size();	}

		int getPeptideLog2RatioSize()
		{ return vdPeptideLog2Ratio.size();	}

	private:

		// utility function for caculating accumulative probility for normal distribution
		double pnorm( double dMean, double dStandardDeviation,
				double dRandomVariable );

		// utility function for calculating log2
		double logBase2( double dNum )
		{ return (double) ( log( dNum ) / log( BASE ) ); }

		// ???
		// input variables
		vector< double > vdPeptideLog2Ratio;
		vector< double > vdPeptideLog2SN;
		
		// Result variables
		double dProteinLog2Ratio;
		double dLowerLimitCI;
		double dUpperLimitCI;
		vector<double> vdLog2RatioBin;
		vector<double> vdLnLikelihood;
		double dMaxLnLikelihood;
		double dLnLikelihoodCutoff;	
};

#endif // PROTEINRATIO_H
