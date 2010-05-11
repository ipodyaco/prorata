
#include "proteinRatio.h"

ProteinRatio::ProteinRatio()
{
	dProteinLog2Ratio = 0;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;
	dMaxLnLikelihood = 0;
	dLnLikelihoodCutoff = 0;
}

ProteinRatio::~ProteinRatio()
{

}

bool ProteinRatio::process( const vector<double> & vdPeptideLog2RatioInput, 
		const vector<double> & vdPeptideLog2SNInput )
{
	int i;
	vdPeptideLog2SN.clear() ;
	vdPeptideLog2Ratio.clear();
	for(  i = 0 ; i < vdPeptideLog2RatioInput.size(); ++i )
	{
		vdPeptideLog2Ratio.push_back( vdPeptideLog2RatioInput[i] );
		vdPeptideLog2SN.push_back( vdPeptideLog2SNInput[i] );
	}

	double dMinLog2SNR = ProRataConfig::getMinLog2SNR();
	double dMaxLog2SNR = ProRataConfig::getMaxLog2SNR();

	double dMinProteinLog2Ratio = ProRataConfig::getMLEMinLog2Ratio();
	double dMaxProteinLog2Ratio = ProRataConfig::getMLEMaxLog2Ratio();

	if ( dMinProteinLog2Ratio > dMaxProteinLog2Ratio )
	{
		cout << "Low and High Ratio values are not appropriate, Please check the session configuration file." << endl;
		return false;
	}

	// if there is only two peptides and their log2 abundance ratio difference is greater than 4
	// then filter them out by setting CI 
	if( vdPeptideLog2Ratio.size() == 2 )
	{
		if( fabs(vdPeptideLog2Ratio[0] - vdPeptideLog2Ratio[1]) > 4 )
		{
			
			dLowerLimitCI = dMinProteinLog2Ratio;
			dUpperLimitCI = dMaxProteinLog2Ratio;
			dProteinLog2Ratio = 0;
			return true;
		}
	}

	double dDiscretizationInterval = ProRataConfig::getLog2RatioDiscretization();
	
	double dSDSlope = ProRataConfig::getSDSlope();
	double dSDIntercept = ProRataConfig::getSDIntercept();
	
	double dMeanSlope = ProRataConfig::getMeanSlope();
	double dMeanIntercept = ProRataConfig::getMeanIntercept();

	if( (dSDIntercept + dSDSlope * dMaxLog2SNR) <= 0 )
	{
		cout << "The predicted standard deviation of the peptide abundance ratios could be negative!" << endl;
		cout << "Change the config <STANDARD_DEVIATION><SLOPE>, <STANDARD_DEVIATION><INTERCEPT> and/or <MAX_LOG2_SNR>!" << endl;
		return false;
	}

	
	double dSmoothingProbabilitySpace = ProRataConfig::getSmoothingProbilitySpace();
	double dLnLikelihoodCutoffOffset = ProRataConfig::getLnLikelihoodCutoffOffset();
	
	/*
	 * set up the array of ratios whose likelihoods are going to be caculated
	 */

	double dCurrentLogRatio = dMinProteinLog2Ratio;

	while( dCurrentLogRatio <= dMaxProteinLog2Ratio)
	{
		vdLog2RatioBin.push_back( dCurrentLogRatio );
		dCurrentLogRatio = dCurrentLogRatio + dDiscretizationInterval;
		// make it exact zero
		if( fabs( dCurrentLogRatio ) < 0.000000001 )
			dCurrentLogRatio = 0;
	}	

	int iRatioBinCount = vdLog2RatioBin.size();

	double dSmoothingProbability = dSmoothingProbabilitySpace / (double) iRatioBinCount;
	
	double dCurrentLnLikelihood = 0.0;
	double dExpectedLogRatio = 0.0;
	double dMean = 0.0;
	double dStandardDeviation = 0.0;
	double dNormProbability = 0.0;
	double dProbability = 0.0;

	for( int n = 0; n < iRatioBinCount; ++n )
	{
		dCurrentLnLikelihood = 0;
		dExpectedLogRatio = vdLog2RatioBin[n];
		for( i = 0; i < vdPeptideLog2Ratio.size(); ++i )
		{
			if(vdPeptideLog2SN[i] < dMinLog2SNR )
				continue;

			// calculate the mean
			dMean=0;
			if( ( ( dMeanSlope * vdPeptideLog2SN[i] + dMeanIntercept ) <= fabs( dExpectedLogRatio ) ) && ( vdPeptideLog2SN[i] < dMaxLog2SNR ) )
			{
				dMean = ( dMeanSlope * vdPeptideLog2SN[i] + dMeanIntercept ) * ( dExpectedLogRatio / fabs( dExpectedLogRatio ) );
			}
			else
			{
				dMean = dExpectedLogRatio;
			}

			// calculate the standard deviation
			dStandardDeviation = 0;
			if(vdPeptideLog2SN[i] < dMaxLog2SNR)
			{
				dStandardDeviation = dSDIntercept + dSDSlope * vdPeptideLog2SN[i];
			}
			else
			{
				dStandardDeviation = dSDIntercept + dSDSlope * dMaxLog2SNR;
			}

			// calculate the probability
			double dRandVar1 = (vdPeptideLog2Ratio.at(i) + (dDiscretizationInterval / 2.0) );
			double dRandVar2 = (vdPeptideLog2Ratio.at(i) - (dDiscretizationInterval / 2.0) );

			dNormProbability = pnorm(dMean, dStandardDeviation, dRandVar1) - pnorm(dMean, dStandardDeviation, dRandVar2);

			dProbability = dNormProbability * ( 1 - dSmoothingProbabilitySpace ) + dSmoothingProbability;

			dCurrentLnLikelihood = dCurrentLnLikelihood + log(dProbability);
				
		}
		vdLnLikelihood.push_back(dCurrentLnLikelihood);
	}

/*
	//Printing stuff
	//cout << "Size of vdLog2RatioBin = " << vdLog2RatioBin.size() << endl;
	//cout << "Size of vdLnLikelihood = " << vdLnLikelihood.size() << endl;

	for( i = 0; i < vdLog2RatioBin.size(); i++ )
	{
		
		cout << "i = " << i << "; vdLnLikelihood = " << vdLnLikelihood.at(i)
			<< "  & vdLog2RatioBin = " << vdLog2RatioBin.at(i) << endl;
			
		//cout << vdLnLikelihood.at(i) << endl;
	}

*/


	dMaxLnLikelihood = ( *(max_element( vdLnLikelihood.begin(), vdLnLikelihood.end() ) ) );

	vector<double> vdSubRatioArray;

	// there could be multiple ratios that have maximum likelihood
	for( i = 0; i < vdLog2RatioBin.size(); i++ )
	{
		if ( vdLnLikelihood.at( i ) == dMaxLnLikelihood )
		{
			vdSubRatioArray.push_back( vdLog2RatioBin.at( i ) );
		}
	}

	
	dProteinLog2Ratio = vdSubRatioArray[0];
	for ( i = 1; i < vdSubRatioArray.size(); i++ )
	{
		if ( fabs( dProteinLog2Ratio ) > fabs( vdSubRatioArray.at( i ) ) )
		{
				dProteinLog2Ratio = vdSubRatioArray.at( i );
		}
	}

//	dLnLikelihoodCutoff = dMaxLnLikelihood - LOGLIKEOFFSET;
	dLnLikelihoodCutoff = dMaxLnLikelihood - dLnLikelihoodCutoffOffset;

	vdSubRatioArray.clear();

	for( i = 0; i < vdLog2RatioBin.size(); i++ )
	{
		if ( vdLnLikelihood.at( i ) > dLnLikelihoodCutoff )
		{
			vdSubRatioArray.push_back( vdLog2RatioBin.at( i ) );
		}
	}

	dLowerLimitCI = * min_element( vdSubRatioArray.begin(), vdSubRatioArray.end() );
	dUpperLimitCI = * max_element( vdSubRatioArray.begin(), vdSubRatioArray.end() );
	
	dLowerLimitCI = dLowerLimitCI - dDiscretizationInterval;
	dUpperLimitCI = dUpperLimitCI + dDiscretizationInterval;

	if( fabs( dLowerLimitCI ) < 0.000000001 )
		dLowerLimitCI = 0;
	if( dLowerLimitCI < dMinProteinLog2Ratio )
	       dLowerLimitCI = dMinProteinLog2Ratio;	

	if( fabs( dUpperLimitCI ) < 0.000000001 )
		dUpperLimitCI = 0;
	if( dUpperLimitCI > dMaxProteinLog2Ratio )
		dUpperLimitCI = dMaxProteinLog2Ratio;

	/*
	sort( vdSubRatioArray.begin(), vdSubRatioArray.end() );

	if ( vdSubRatioArray.size() > 0 )
	{
		dLowerLimitCI = vdSubRatioArray.at( 0 );
		dUpperLimitCI = vdSubRatioArray.at( vdSubRatioArray.size() - 1 );
	}
	*/

	return true;
}

double ProteinRatio::pnorm( double dMean, double dStandardDeviation,
				double dRandomVariable )
{
	double dZScore = ( dRandomVariable - dMean ) / dStandardDeviation ;

	double dProbability;

	/* 
	 * this implementation uses a Linux-specific function erfc
	 */
	 dProbability = PNORM_COEFF1 * erfc( -dZScore / sqrt( PNORM_COEFF2 ) );

	/*
	 * this implementation uses a simple approximation
	 */
/*
	  if( dZScore > 6.0 ) { return 1.0; }
	  if( dZScore < -6.0 ) { return 0.0; }

	  double b1 = 0.31938153;
	  double b2 = -0.356563782;
	  double b3 = 1.781477937;
	  double b4 = -1.821255978;
	  double b5 = 1.330274429;
	  double p = 0.2316419;
	  double c2 = 0.3989423;

	  double a = fabs(dZScore);
	  double t = 1.0/(1.0+a*p);
	  double b = c2*exp((-dZScore)*(dZScore/2.0));
	  dProbability = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	  dProbability = 1.0 - b*dProbability;
	  if( dZScore < 0.0 ) dProbability = 1.0 - dProbability;
*/
	  
	  return dProbability;
}


