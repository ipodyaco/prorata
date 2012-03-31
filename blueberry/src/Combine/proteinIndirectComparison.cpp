
#include "proteinIndirectComparison.h"

double ProteinIndirectComparison::dMaxCIwidth = 3;
double ProteinIndirectComparison::dLnLikehoodCutOffset = 1.96;

ProteinIndirectComparison::ProteinIndirectComparison()
{
	bValidity = false;
	sLocus = "";
	dLog2Ratio = 0;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;
	sName = "";
	sNumerator = "";
	sDenominator = "";
	sDirectComparisonName0 = "";
	sDirectComparisonName1 = "";
	dMLEMinLog2Ratio = 0;
	dMLEMaxLog2Ratio = 0;
	dLog2RatioDiscretization = 0.1;
}

ProteinIndirectComparison::~ProteinIndirectComparison()
{
	// destructor
}

bool ProteinIndirectComparison::testDirectComparison(
		ProteinDirectComparison directComparison0, ProteinDirectComparison directComparison1)
{
	if( !calCombinationType( directComparison0.getNumerator(), directComparison0.getDenominator(),
				directComparison1.getNumerator(), directComparison1.getDenominator() ) )
	{
		cout << "ERROR: cannot properly match the numerator and denominator of the two direct comparisons and the indirect comparison." << endl;
		return false;
	}

	if( directComparison0.getMLEMinLog2Ratio() != directComparison1.getMLEMinLog2Ratio() )
	{
		cout << "ERROR: the two direct comparison must have the same minimum protein log2ratio." << endl;
		return false;
	}
	
	if( directComparison0.getMLEMaxLog2Ratio() != directComparison1.getMLEMaxLog2Ratio() )
	{
		cout << "ERROR: the two direct comparison must have the same maximum protein log2ratio." << endl;
		return false;
	}
	
	if( directComparison0.getLog2RatioDiscretization() != directComparison1.getLog2RatioDiscretization() )
	{
		cout << "ERROR: the two direct comparison must have the same protein Log2RatioDiscretization." << endl;
		return false;
	}

	return true;

}

bool ProteinIndirectComparison::runQuantification( string sInputLocus,
		ProteinDirectComparison directComparison0, ProteinDirectComparison directComparison1)
{

	bValidity = false;
	dLog2Ratio = 0;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;
	if( !( directComparison0.isValid() && directComparison1.isValid() ) )
	{
		// both direct comparisons for this locus have to be valid
		return true;
	}

	
	if( !calCombinationType( directComparison0.getNumerator(), directComparison0.getDenominator(),
				directComparison1.getNumerator(), directComparison1.getDenominator() ) )
	{
		cout << "ERROR: cannot properly match the numerator and denominator of the two direct comparisons and the indirect comparison." << endl;
		return false;
	}
	
	vector<double> vdLog2RatioDC0;
	vector<double> vdLnLikelihoodDC0;	
	vector<double> vdLog2RatioDC1;
	vector<double> vdLnLikelihoodDC1;

	if( !directComparison0.getProfileLikelihoodCurve(vdLog2RatioDC0, vdLnLikelihoodDC0) )
		return false;
	
	if( !directComparison1.getProfileLikelihoodCurve(vdLog2RatioDC1, vdLnLikelihoodDC1) )
		return false;

	if( vdLog2RatioDC0.size() != vdLog2RatioDC1.size() )
		return false;
	
	int i;
	for( i = 0; i < vdLog2RatioDC0.size(); ++i )
	{
		if( fabs(vdLog2RatioDC0[i] - vdLog2RatioDC1[i]) > 0.001 )
			return false;
	}

	
	// combine every log2ratio from comparison 0 with every log2ratio from comparison 1 
	vector<double> vdLog2RatioAll;
	vector<double> vdLnLikelihoodAll;	
	int j;
	for( i = 0; i < vdLog2RatioDC0.size(); ++i )
	{
		for( j = 0; j < vdLog2RatioDC1.size(); ++j )
		{
			vdLog2RatioAll.push_back( combineLog2Ratio(vdLog2RatioDC0[i], vdLog2RatioDC1[j]) );
			vdLnLikelihoodAll.push_back( (vdLnLikelihoodDC0[i] + vdLnLikelihoodDC1[j]) );
		}
	}

	// collapse the log2ratio into a non-redundant set of log2ratio
	double dCurrentLog2Ratio = 0;
	double dCurrentLnLikelihood = 0;
	double dOverMaxLog2Ratio = 10000000000.0;
	bool bFirstLog2Ratio = true;
	vdLog2Ratio.clear();
	vdLnLikelihoodAll.clear();

	dCurrentLog2Ratio = ( *(min_element( vdLog2RatioAll.begin(), vdLog2RatioAll.end() ) ) );
	while( dCurrentLog2Ratio < dOverMaxLog2Ratio )
	{
		bFirstLog2Ratio = true;
		dCurrentLnLikelihood = 0;
	//	cout << "dCurrentLog2Ratio = "  << dCurrentLog2Ratio << endl;
		for( i = 0; i < vdLog2RatioAll.size(); ++i )
		{
			if( fabs(vdLog2RatioAll[i] - dCurrentLog2Ratio) < 0.00001 )
			{
				// compute ln likelihood
				if(bFirstLog2Ratio)
				{
					dCurrentLnLikelihood = vdLnLikelihoodAll[i];
					bFirstLog2Ratio = false;
				}
				else
				{
					// Wrong implementation: exp(vdLnLikelihoodAll[i]) can easily be out of range for double
					// dCurrentLnLikelihood = log( exp(dCurrentLnLikelihood) + exp(vdLnLikelihoodAll[i]) );
					// if y > x, log(exp(x)+exp(y)) = y + log(exp(x-y)+1)
					if( dCurrentLnLikelihood > vdLnLikelihoodAll[i])
					{
						dCurrentLnLikelihood = dCurrentLnLikelihood + log( 1 + exp(vdLnLikelihoodAll[i]-dCurrentLnLikelihood) );
					}
					else
					{
						dCurrentLnLikelihood = vdLnLikelihoodAll[i] + log( 1 + exp(dCurrentLnLikelihood-vdLnLikelihoodAll[i]) );
					}
					
				}
				vdLog2RatioAll[i] = dOverMaxLog2Ratio;
				//cout << "vdLog2RatioAll[i] = "  << vdLog2RatioAll[i] << endl;
			}
		}
		vdLog2Ratio.push_back( dCurrentLog2Ratio );
		vdLnLikelihood.push_back( dCurrentLnLikelihood );
		dCurrentLog2Ratio = ( *(min_element( vdLog2RatioAll.begin(), vdLog2RatioAll.end() ) ) );
	}			

	//cout << "vdLog2Ratio.size() = " << vdLog2Ratio.size() << endl;
	//for( int x = 0; x < vdLog2Ratio.size(); ++x )
	//{
	//	cout << "##  " << vdLog2Ratio[x] << " === " << vdLnLikelihood[x] << endl; 
	//}
	
	vdLog2RatioAll.clear();
	vdLnLikelihoodAll.clear();

	/*
	 * calculate the maximum likelihood estimate and confidence interval from profile likelihood curve
	 * this is copied from proteinDirectComparison.cpp
	 * changed variable names:  dMLEMinLog2Ratio dMLEMaxLog2Ratio dLog2RatioDiscretization
	 */
	
	dMLEMinLog2Ratio = ( *(min_element( vdLog2Ratio.begin(), vdLog2Ratio.end() ) ) );
	dMLEMaxLog2Ratio = ( *(max_element( vdLog2Ratio.begin(), vdLog2Ratio.end() ) ) );
	dLog2RatioDiscretization = directComparison0.getLog2RatioDiscretization();

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

	/*
	 * Determine the validity of this protein according to the filterin criteria
	 */
	
	// CI width filtering
	double dWidthCI = dUpperLimitCI - dLowerLimitCI;
	if( dWidthCI < dMaxCIwidth )
		bValidity = true;
	else
		bValidity = false;

	

	return true;
}


void ProteinIndirectComparison::eraseProfileLikelihoodCurve()
{
	vdLog2Ratio.clear();
	vdLnLikelihood.clear();
}

void ProteinIndirectComparison::setName( string sNameConfig )
{
	sName = sNameConfig;
}

void ProteinIndirectComparison::setNumerator( string sNumeratorConfig )
{
	sNumerator = sNumeratorConfig;
}

void ProteinIndirectComparison::setDenominator( string sDenominatorConfig )
{
	sDenominator = sDenominatorConfig;
}

void ProteinIndirectComparison::setDirectComparisonName0( string sDirectComparisonName0Config)
{
	sDirectComparisonName0 = sDirectComparisonName0Config;
}

void ProteinIndirectComparison::setDirectComparisonName1( string sDirectComparisonName1Config )
{
	sDirectComparisonName1 = sDirectComparisonName1Config;
}

string ProteinIndirectComparison::getName()
{
	return sName;
}

string ProteinIndirectComparison::getNumerator()
{
	return sNumerator;
}

string ProteinIndirectComparison::getDenominator()	
{
	return sDenominator;
}

string ProteinIndirectComparison::getDirectComparisonName0()
{
	return sDirectComparisonName0;
}

string ProteinIndirectComparison::getDirectComparisonName1()
{
	return sDirectComparisonName1;
}	

void ProteinIndirectComparison::getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef )
{
	bValidityRef = bValidity;
	dLog2RatioRef = dLog2Ratio;
	dLowerLimitCIRef = dLowerLimitCI;
	dUpperLimitCIRef = dUpperLimitCI;
}
		
bool ProteinIndirectComparison::isValid()
{
	return bValidity;
}	

bool ProteinIndirectComparison::calCombinationType(string sDC0Numerator, string sDC0Denominator,
		string sDC1Numerator, string sDC1Denominator)
{
	if( ( (sNumerator == sDC0Numerator) && (sDenominator == sDC1Numerator) && (sDC0Denominator == sDC1Denominator) )
		|| ( (sNumerator == sDC1Denominator) && (sDenominator == sDC0Denominator) && (sDC0Numerator == sDC1Numerator) ) )
	{
		iCombinationType = 1;
		return true;
	}
	else if( ( (sNumerator == sDC0Denominator) && (sDenominator == sDC1Denominator) && (sDC0Numerator == sDC1Numerator) )
		|| ( (sNumerator == sDC1Numerator) && (sDenominator == sDC0Numerator) && (sDC0Denominator == sDC1Denominator) ) )
	{
		iCombinationType = 2;
		return true;
	}
	else if( ( (sNumerator == sDC0Numerator) && (sDenominator == sDC1Denominator) && (sDC0Denominator == sDC1Numerator) )
		|| ( (sNumerator == sDC1Numerator) && (sDenominator == sDC0Denominator) && (sDC0Numerator == sDC1Denominator) ) )
	{
		iCombinationType = 3;
		return true;
	}
	else if( ( (sNumerator == sDC0Denominator) && (sDenominator == sDC1Numerator) && (sDC0Numerator == sDC1Denominator) )
		|| ( (sNumerator == sDC1Denominator) && (sDenominator == sDC0Numerator) && (sDC0Denominator == sDC1Numerator) ) )
	{
		iCombinationType = 4;
		return true;
	}
	else
	{
		iCombinationType = 0;
		return false;
	}

}

double ProteinIndirectComparison::combineLog2Ratio( double dLog2RatioDC0, double dLog2RatioDC1 )
{
	if( iCombinationType == 1 )
		return (dLog2RatioDC0 - dLog2RatioDC1);
	else if( iCombinationType == 2 )
		return (-(dLog2RatioDC0 - dLog2RatioDC1));
	else if( iCombinationType == 3 )
		return (dLog2RatioDC0 + dLog2RatioDC1);
	else if( iCombinationType == 4 )
		return (-(dLog2RatioDC0 + dLog2RatioDC1));
	else
		return 0;
}




