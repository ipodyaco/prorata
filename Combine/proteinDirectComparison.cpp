
#include "proteinDirectComparison.h"

int ProteinDirectComparison::iValidReplicateCICutoff = 2;
int ProteinDirectComparison::iMinimumPeptideNumber = 2;
double ProteinDirectComparison::dMaxCIwidth = 3.0;
double ProteinDirectComparison::dLnLikehoodCutOffset = 1.96;
bool ProteinDirectComparison::bRequireCIoverlapped = true;
double ProteinDirectComparison::dMaxLogRatioDifference = 2.0;

ProteinDirectComparison::ProteinDirectComparison()
{
	sName = "";
	sNumerator = "";
	sDenominator = "";
	bValidity = false;
	sLocus = "";
	dLog2Ratio = 0;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;
	iValidReplicateNumber = 0;
	iValidPeptideNumber = 0;
	iQuantifiedPeptides = 0;
	dMLEMinLog2Ratio = 0;
	dMLEMaxLog2Ratio = 0;
	dLog2RatioDiscretization = 0.1;
}

ProteinDirectComparison::~ProteinDirectComparison()
{
	// destructor
}

bool ProteinDirectComparison::runQuantification( string sInputLocus )
{
	sLocus = sInputLocus;
	dLog2Ratio = 0;
	bValidity = false;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;

	iValidReplicateNumber = 0;
	iValidPeptideNumber = 0;

	/*
	 *  run quantification for all replicates
	 */
	int i = 0;
	for( i = 0; i < vReplicate.size(); ++i )
	{
		if( !vReplicate[i].runQuantification( sLocus ) )
		{
			cout << "WARNNING: cannot run quantification for locus " << sLocus << " in Replicate " << vReplicate[i].getName() << endl;
		}
	}
	

	// setup vdLog2Ratio and vdLnLikelihood
	vdLog2Ratio.clear();
	vdLnLikelihood.clear();
	double dCurrentLogRatio = dMLEMinLog2Ratio;
	while( dCurrentLogRatio <= dMLEMaxLog2Ratio)
	{
		vdLog2Ratio.push_back( dCurrentLogRatio );
		vdLnLikelihood.push_back(0);
		dCurrentLogRatio = dCurrentLogRatio + dLog2RatioDiscretization;
		// make it exact zero
		if( fabs( dCurrentLogRatio ) < 0.000000001 )
			dCurrentLogRatio = 0;
	}

	if(vdLog2Ratio.size() < 3)
	{
		cout << "ERROR: Problem with maximum log2ratio, minimum log2ratio or discretization interval at Direct Comparison." << endl;
		bValidity = false;
		return false;

	}
	

	/*
	 * combine the profile likelihood curves from the replicates
	 */

	vector<double> vdLog2RatioTemp;
	vector<double> vdLnLikelihoodTemp;
	

	
	bool bFlagFirstValidReplicate = true;
	int j;
	unsigned int k;
	double dTempLog2RatioDifference = 1000000.0;
	double dTempLnLikelihood = 0;
	for( i = 0; i < vReplicate.size(); ++i )
	{
		if( vReplicate[i].isValid() )
		{
			bFlagFirstValidReplicate = false;
			vReplicate[i].getProfileLikelihoodCurve( vdLog2RatioTemp, vdLnLikelihoodTemp );

			for( j = 0; j < vdLog2Ratio.size(); ++j )
			{
				dTempLog2RatioDifference = 1000000.0;
				for( k = 0; k < vdLog2RatioTemp.size(); ++k )
				{
					if( dTempLog2RatioDifference > fabs(vdLog2Ratio[j] - vdLog2RatioTemp[k]))
					{
						dTempLog2RatioDifference =  fabs(vdLog2Ratio[j] - vdLog2RatioTemp[k]);
						dTempLnLikelihood = vdLnLikelihoodTemp[k];
					}
				}
				vdLnLikelihood[j] = vdLnLikelihood[j] + dTempLnLikelihood;
			}

		}
	}
	/*
	for( i = 0; i < vdLnLikelihood.size(); i++ )
	{
		cout << vdLog2Ratio[i] << " - " << vdLnLikelihood[i] << endl;

	}
	*/
	if(bFlagFirstValidReplicate)
	{
		// this locus is not found in any replicates
		bValidity = false;
		return true;
	}

	// calculate iQuantifiedPeptides
	iQuantifiedPeptides = 0;
	for( i = 0; i < vReplicate.size(); ++i )
	{
		iQuantifiedPeptides = iQuantifiedPeptides + vReplicate[i].getQuantifiedPeptides();
	}

	/*
	 * calculate the maximum likelihood estimate and confidence interval from profile likelihood curve
	 * this is copied from proteinRatio.cpp
	 * changed variable names: dLog2Ratio dMLEMinLog2Ratio dMLEMaxLog2Ratio dLog2RatioDiscretization
	 */
	double dMaxLnLikelihood = ( *(max_element( vdLnLikelihood.begin(), vdLnLikelihood.end() ) ) );
	vector<double> vdSubRatioArray;

	// there could be multiple ratios that have maximum likelihood
	for( i = 0; i < vdLog2Ratio.size(); i++ )
	{
		if ( vdLnLikelihood.at( i ) == dMaxLnLikelihood )
		{
			vdSubRatioArray.push_back( vdLog2Ratio.at( i ) );
		}
	}

	dLog2Ratio = vdSubRatioArray[0];
	for ( i = 1; i < vdSubRatioArray.size(); i++ )
	{
		if ( fabs( dLog2Ratio ) > fabs( vdSubRatioArray.at( i ) ) )
		{
				dLog2Ratio = vdSubRatioArray.at( i );
		}
	}

	double dLnLikelihoodCutoff = dMaxLnLikelihood - dLnLikehoodCutOffset;

	vdSubRatioArray.clear();

	for( i = 0; i < vdLog2Ratio.size(); i++ )
	{
		if ( vdLnLikelihood.at( i ) > dLnLikelihoodCutoff )
		{
			vdSubRatioArray.push_back( vdLog2Ratio.at( i ) );
		}
	}

	dLowerLimitCI = * min_element( vdSubRatioArray.begin(), vdSubRatioArray.end() );
	dUpperLimitCI = * max_element( vdSubRatioArray.begin(), vdSubRatioArray.end() );
	
	dLowerLimitCI = dLowerLimitCI - dLog2RatioDiscretization;
	dUpperLimitCI = dUpperLimitCI + dLog2RatioDiscretization;

	if( fabs( dLowerLimitCI ) < 0.000000001 )
		dLowerLimitCI = 0;
	if( dLowerLimitCI < dMLEMinLog2Ratio )
	       dLowerLimitCI = dMLEMinLog2Ratio;	

	if( fabs( dUpperLimitCI ) < 0.000000001 )
		dUpperLimitCI = 0;
	if( dUpperLimitCI > dMLEMaxLog2Ratio )
		dUpperLimitCI = dMLEMaxLog2Ratio;

//	cout << "Direct Comparison dLog2Ratio = " << dLog2Ratio << " dUpperLimitCI = " << dUpperLimitCI << " dLowerLimitCI " << dLowerLimitCI << endl;


	/*
	 * Determine the validity of this protein according to the filterin criteria
	 */
	
	// CI width filtering
	double dWidthCI = dUpperLimitCI - dLowerLimitCI;
	if( dWidthCI > dMaxCIwidth )
	{
		bValidity = false;
		return true;
	}

	j = 0;
	bool bValidity0 = false;
	double dLogRatio0 = 0;
	double dCI0upper = 0;
	double dCI0lower = 0;
	bool bValidity1 = false;
	double dLogRatio1 = 0;
	double dCI1upper = 0;
	double dCI1lower = 0;
	bool bValidReplicate = false;
	iValidPeptideNumber = 0;
	iValidReplicateNumber = 0;
	for( i = 0; i < vReplicate.size(); ++i )
	{
		if(!vReplicate[i].isValid())
			continue;

		// Implementation 3: no peptide number requirement
		/*
		// add up valid peptides in each valid replicate
		iValidPeptideNumber += vReplicate[i].getPeptideNumberWithinCI( dLowerLimitCI, dUpperLimitCI );
		*/

		// determine valid replicate number
		vReplicate[i].getEstimates(bValidity0, dLogRatio0, dCI0lower, dCI0upper);


		//Implementation 3: CI overlap or low log2ratio error
		if( bIsCIOverlapped(dCI0lower, dCI0upper, dLowerLimitCI, dUpperLimitCI) || fabs(dLogRatio0-dLog2Ratio) < dMaxLogRatioDifference )
		{
			iValidReplicateNumber++;
		}


		// Implementation 2 : require the CI overlap between a replicate and the combined CI
		/*
		if( bIsCIOverlapped(dCI0lower, dCI0upper, dLowerLimitCI, dUpperLimitCI) && fabs(dLogRatio0-dLog2Ratio) < dMaxLogRatioDifference )
		{
			iValidReplicateNumber++;
		}
		*/

		//Implementation 1 : require the CI overlap between replicates
		/*
		if( !bIsCIOverlapped(dCI0lower, dCI0upper, dLowerLimitCI, dUpperLimitCI) )
		{
			continue;
		}

		bValidReplicate = false;
		for( j = 0; j < vReplicate.size(); ++j )
		{
			if( j != i && vReplicate[j].isValid() )
			{
				vReplicate[j].getEstimates(bValidity1, dLogRatio1, dCI1lower, dCI1upper);
				if( bRequireCIoverlapped )
				{
					if( bIsCIOverlapped(dCI0lower, dCI0upper, dCI1lower, dCI1upper)
							&& ( fabs( dLogRatio0 - dLogRatio1 ) < dMaxLogRatioDifference))
					{
						bValidReplicate = true;
						break;
					}
				}
				else
				{
					if( fabs(dLogRatio0-dLogRatio1) < dMaxLogRatioDifference)
					{
						bValidReplicate = true;
						break;
					}
				}
			}
		}
		if( bValidReplicate )
		{
			iValidReplicateNumber++;
		}
		*/

	}
		
	
	//if( iValidPeptideNumber < iMinimumPeptideNumber || iValidReplicateNumber < iValidReplicateCICutoff  )
	if( iValidReplicateNumber < iValidReplicateCICutoff  )
	{
		bValidity = false;
	}
	else
	{
		bValidity = true;
	}

	
	return true;
}

void ProteinDirectComparison::eraseProfileLikelihoodCurve()
{
	vdLog2Ratio.clear();
	vdLnLikelihood.clear();
	for( int i = 0; i < vReplicate.size(); ++i )
	{
		vReplicate[i].eraseProfileLikelihoodCurve();
	}
}

bool ProteinDirectComparison::getProfileLikelihoodCurve(vector<double> & vdLog2RatioRef, vector<double> & vdLnLikelihoodRef )
{
	vdLog2RatioRef.clear();
	vdLnLikelihoodRef.clear();
	
	if( vdLog2Ratio.size() == 0 || vdLnLikelihood.size() == 0 )
	{
		return false;
	}
	
	vdLog2RatioRef = vdLog2Ratio;
	vdLnLikelihoodRef = vdLnLikelihood;
	return true;
}

const vector< ProteinReplicate > & ProteinDirectComparison::getProteinReplicateVector()
{
	return vReplicate;
}

void ProteinDirectComparison::setName( string sNameConfig )
{
	sName = sNameConfig;
}

void ProteinDirectComparison::setNumerator( string sNumeratorConfig )
{
	sNumerator = sNumeratorConfig;
}

void ProteinDirectComparison::setDenominator( string sDenominatorConfig )
{
	sDenominator = sDenominatorConfig;
}

string ProteinDirectComparison::getName()
{
	return sName;
}

string ProteinDirectComparison::getNumerator()
{
	return sNumerator;
}

string ProteinDirectComparison::getDenominator()
{
	return sDenominator;
}
/*
double  ProteinDirectComparison::getMLEMinLog2Ratio()
{
	return dMLEMinLog2Ratio;
}

double  ProteinDirectComparison::getMLEMaxLog2Ratio()
{
	return dMLEMaxLog2Ratio;
}
*/
double  ProteinDirectComparison::getLog2RatioDiscretization()
{
	return dLog2RatioDiscretization;
}

bool ProteinDirectComparison::addReplicate( ProteinReplicate tempProteinReplicate )
{
	double dMLEMinLog2RatioTemp = 0;
	double dMLEMaxLog2RatioTemp = 0;
	double dLog2RatioDiscretizationTemp = 0.1;

	if( !tempProteinReplicate.checkRatio(sNumerator, sDenominator) )
	{
		cout << "ERROR: Inconsistant numerator and denominator between Direct Comparison " << sName << " and Replicate " << tempProteinReplicate.getName() << endl;
		return false;
	}

	tempProteinReplicate.computeLog2RatioSetting( dMLEMinLog2RatioTemp, dMLEMaxLog2RatioTemp, dLog2RatioDiscretizationTemp );

	if( vReplicate.size() == 0 )
	{
		// the first replicate
		dLog2RatioDiscretization = dLog2RatioDiscretizationTemp;
		dMLEMinLog2Ratio = dMLEMinLog2RatioTemp;
		dMLEMaxLog2Ratio = dMLEMaxLog2RatioTemp;		
	}
	else
	{
		if(dLog2RatioDiscretization != dLog2RatioDiscretizationTemp)
		{
			cout << "ERROR: all replicates to be combined must have the same <LOG2_RATIO_DISCRETIZATION> in their ProRataConfig.xml! " << endl;
			return false;
		}
		if( dMLEMinLog2Ratio > dMLEMinLog2RatioTemp )
		{
			dMLEMinLog2Ratio = dMLEMinLog2RatioTemp;
		}
		if( dMLEMaxLog2Ratio < dMLEMaxLog2RatioTemp )
		{
			dMLEMaxLog2Ratio = dMLEMaxLog2RatioTemp;
		}		
	}
	vReplicate.push_back( tempProteinReplicate );
	return true;
}

bool ProteinDirectComparison::bIsCIOverlapped(double dCI0lower, double dCI0upper, double dCI1lower, double dCI1upper)
{
	if( (dCI0lower <= dCI1upper) && (dCI1lower <= dCI0upper) )
		return true;
	else
		return false;
}


string ProteinDirectComparison::getLocus()
{
	return sLocus;
}

bool ProteinDirectComparison::isValid()
{
	return bValidity;
}

void ProteinDirectComparison::getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef )
{
	bValidityRef = bValidity;
	dLog2RatioRef = dLog2Ratio;
	dLowerLimitCIRef = dLowerLimitCI;
	dUpperLimitCIRef = dUpperLimitCI;
}

int ProteinDirectComparison::getValidReplicateNumber()
{
	return iValidReplicateNumber;
}

int ProteinDirectComparison::getValidPeptideNumber()
{
	return iValidPeptideNumber;
}

int ProteinDirectComparison::getQuantifiedPeptides()
{
	return iQuantifiedPeptides;
}


