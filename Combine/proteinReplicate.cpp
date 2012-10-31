#include "proteinReplicate.h"

ProteinReplicate::ProteinReplicate()
{
	sName = "";
	sNormalizationMethod = "";
	dNormalizationValue = 0;	
	pProteomeInfo = NULL;

	bRemoveAmbiguousPeptide = true;
	iMinPeptideNumber = 0;
	dMaxCIwidth = 0;
	dMinLog2SNR = 0;
	dMaxLog2SNR = 0;
	dMLEMinLog2Ratio = 0;
	dMLEMaxLog2Ratio = 0;
	dLog2RatioDiscretization = 0.1;
	dSDSlope = 0;
	dSDIntercept = 0;
	dMeanSlope = 0;
	dMeanIntercept = 0; 
	dSmoothingProbSpace = 0;
	dLnLikelihoodCutoffOffset = 0;

	sOriginalNumerator = "";
	sOriginalDenominator = "";
	bReverseRatio = false;
	
	sLocus = "";
	pProteinInfo = NULL;

	bValidity = false;
	iQuantifiedPeptides = 0;

	dLog2Ratio = 0;
	dLowerLimitCI = 0;
	dUpperLimitCI = 0;

}

ProteinReplicate::~ProteinReplicate()
{
	// destructor
}

void ProteinReplicate::setName( string sNameConfig )
{
	sName = sNameConfig;
}

void ProteinReplicate::setProteomeInfo( ProteomeInfo * pProteomeInfoConfig , string sNormalizationMethodConfig, double dNormalizationValueConfig )
{
	int iDiscretizationUnits = 0;


	pProteomeInfo = pProteomeInfoConfig;

	bRemoveAmbiguousPeptide	= ProRataConfig::getRemoveAmbiguousPeptides();		
	iMinPeptideNumber	= ProRataConfig::getMinPeptideNumber();			
	dMaxCIwidth		= ProRataConfig::getMaxCIwidth();				
	dMinLog2SNR		= ProRataConfig::getMinLog2SNR();				
	dMaxLog2SNR		= ProRataConfig::getMaxLog2SNR();				
	dMLEMinLog2Ratio	= ProRataConfig::getMLEMinLog2Ratio();			
	dMLEMaxLog2Ratio	= ProRataConfig::getMLEMaxLog2Ratio();			
//	dLog2RatioDiscretization	= ProRataConfig::getLog2RatioDiscretization();

	//  hardcode dLog2RatioDiscretization to 0.1	
	dLog2RatioDiscretization = 0.1;
	
	// round dMLEMinLog2Ratio and dMLEMaxLog2Ratio to units of dLog2RatioDiscretization
	iDiscretizationUnits = (int)( dMLEMinLog2Ratio / dLog2RatioDiscretization + 0.5 );
	dMLEMinLog2Ratio = dLog2RatioDiscretization * iDiscretizationUnits;

	iDiscretizationUnits = (int)( dMLEMaxLog2Ratio / dLog2RatioDiscretization + 0.5 );
	dMLEMaxLog2Ratio = dLog2RatioDiscretization * iDiscretizationUnits;
	
	dSDSlope		= ProRataConfig::getSDSlope();				
	dSDIntercept		= ProRataConfig::getSDIntercept();			
	dMeanSlope		= ProRataConfig::getMeanSlope();				
	dMeanIntercept		= ProRataConfig::getMeanIntercept();			
	dSmoothingProbSpace	= ProRataConfig::getSmoothingProbilitySpace();		
	dLnLikelihoodCutoffOffset	= ProRataConfig::getLnLikelihoodCutoffOffset();
	sOriginalNumerator = ProRataConfig::getNumeratorIsotopologue();
	sOriginalDenominator = ProRataConfig::getDenominatorIsotopologue();	

	sNormalizationMethod = sNormalizationMethodConfig;
	dNormalizationValue = dNormalizationValueConfig;
	// round the dNormalizationValue to units of dLog2RatioDiscretization
	if( dNormalizationValue >= 0 )
		iDiscretizationUnits = (int)( dNormalizationValue / dLog2RatioDiscretization + 0.5 );
	else
		iDiscretizationUnits = (int)( dNormalizationValue / dLog2RatioDiscretization - 0.5 );
	dNormalizationValue = dLog2RatioDiscretization * iDiscretizationUnits;
//	cout << " dNormalizationValue = " << dNormalizationValue << endl;

}

bool ProteinReplicate::checkRatio(string sDirectComparisonNumerator, string sDirectComparisonDenominator)
{
	if(sOriginalNumerator == sDirectComparisonNumerator && 
	   sOriginalDenominator == sDirectComparisonDenominator )
	{
		bReverseRatio = false;
		return true;
	}
	else if( sOriginalDenominator == sDirectComparisonNumerator && 
	         sOriginalNumerator == sDirectComparisonDenominator )
	{
		bReverseRatio = true;
		return true;
	}
	else
	{
		cout << "ERROR: Replicate Numerator : Denominator = " << sOriginalNumerator << " : " << sOriginalDenominator << 
		       	" Direct Comparison Numerator : Denominator = " << sDirectComparisonNumerator << " : " << sDirectComparisonDenominator << endl; 
		return false;
	}	
}

bool ProteinReplicate::runQuantification( string sInputLocus )
{
	// creat the CONFIG environment for this replicate
	ProRataConfig::setRemoveAmbiguousPeptides(bRemoveAmbiguousPeptide);
	ProRataConfig::setMinPeptideNumber(iMinPeptideNumber);
	ProRataConfig::setMaxCIwidth(dMaxCIwidth); 
	ProRataConfig::setMinLog2SNR(dMinLog2SNR); 
	ProRataConfig::setMaxLog2SNR(dMaxLog2SNR); 
	ProRataConfig::setMLEMinLog2Ratio(dMLEMinLog2Ratio);
	ProRataConfig::setMLEMaxLog2Ratio(dMLEMaxLog2Ratio);
	ProRataConfig::setLog2RatioDiscretization(dLog2RatioDiscretization);
	ProRataConfig::setSDSlope(dSDSlope); 
	ProRataConfig::setSDIntercept(dSDIntercept); 
	ProRataConfig::setMeanSlope(dMeanSlope); 
	ProRataConfig::setMeanIntercept(dMeanIntercept); 
	ProRataConfig::setSmoothingProbilitySpace(dSmoothingProbSpace); 
	ProRataConfig::setLnLikelihoodCutoffOffset(dLnLikelihoodCutoffOffset); 

	sLocus = sInputLocus;

	ProteinRatio * pProteinRatio = new ProteinRatio;

	vector< ProteinInfo * > vpProteinInfoSearch;
	vpProteinInfoSearch = pProteomeInfo->getProteinInfo4Locus( sLocus );

	if(vpProteinInfoSearch.size() != 1)
	{
		bValidity = false;
		dLog2Ratio = 0;
		dLowerLimitCI = 0;
		dUpperLimitCI = 0;
		vdLog2Ratio.clear();
		vdLnLikelihood.clear();
		
		// this protein is not quantified in this replicate
		if( vpProteinInfoSearch.size() == 0 )
		{
			// cout << "WARNNING: cannot find locus " << sLocus << " in replicate " << sName << endl;
			return true;
		}

		// this locus occurs more than once, this would be odd
		if( vpProteinInfoSearch.size() > 1 )
		{
			cout << "WARNNING: there are more than one locus named " << sLocus << " in replicate " << sName << endl;
			return true;
		}
	}

	pProteinInfo = vpProteinInfoSearch[0];

	// calculate pProteinRatio	
	if( !pProteomeInfo->getProteinRatio( pProteinInfo, pProteinRatio ) )
	{
		cout << "WARNING: cannot calculate profile likelihood curve for locus " << sLocus << endl;
		return true;
	}

	// retrieve estimation results from pProteinRatio
	vdLog2Ratio = pProteinRatio->getLog2RatioBin();
	vdLnLikelihood = pProteinRatio->getLnLikelihood();
	
	dLog2Ratio = pProteinRatio->getProteinLog2Ratio();
	dLowerLimitCI = pProteinRatio->getLowerLimitCI();
	dUpperLimitCI = pProteinRatio->getUpperLimitCI();

	// normalization
	unsigned int i;
	if( sNormalizationMethod == "Plus" )
	{
		for( i = 0; i < vdLog2Ratio.size(); i++)
		{
			vdLog2Ratio[i] = vdLog2Ratio[i] + dNormalizationValue;
		}

		dLog2Ratio = dLog2Ratio + dNormalizationValue;
		dLowerLimitCI = dLowerLimitCI + dNormalizationValue;
		dUpperLimitCI = dUpperLimitCI + dNormalizationValue;		

	}
	else if ( sNormalizationMethod == "Minus" ) 
	{
		for( i = 0; i < vdLog2Ratio.size(); i++)
		{
			vdLog2Ratio[i] = vdLog2Ratio[i] - dNormalizationValue;
		}

		dLog2Ratio = dLog2Ratio - dNormalizationValue;
		dLowerLimitCI = dLowerLimitCI - dNormalizationValue;
		dUpperLimitCI = dUpperLimitCI - dNormalizationValue;
	}
	else
	{
		// no normalization
	}

	if( fabs( dLowerLimitCI ) < 0.000000001 )
		dLowerLimitCI = 0;
	if( fabs( dUpperLimitCI ) < 0.000000001 )
		dUpperLimitCI = 0;
	if( fabs( dLog2Ratio ) < 0.000000001 )
		dLog2Ratio = 0;

	// reverse the ratio, if needed
	if( bReverseRatio )
	{	
		for( i = 0; i < vdLog2Ratio.size(); i++)
		{
			vdLog2Ratio[i] = -vdLog2Ratio[i];
		}
		reverse( vdLog2Ratio.begin(), vdLog2Ratio.end() );
		reverse( vdLnLikelihood.begin(), vdLnLikelihood.end() );

		dLog2Ratio = -dLog2Ratio;

		double dTemp = dUpperLimitCI;	
		dUpperLimitCI = -dLowerLimitCI ;
		dLowerLimitCI = -dTemp;		
	}

	// validity is retrieved from pProteinInfo read from QPR file
	bValidity = pProteinInfo->getValidity();
	iQuantifiedPeptides = pProteinInfo->getQuantifiedPeptides();
	delete pProteinRatio;

//	cout << " dLog2Ratio = " << dLog2Ratio << " dUpperLimitCI = " << dUpperLimitCI << " dLowerLimitCI " << dLowerLimitCI << endl;
	
	return true;
}

void ProteinReplicate::eraseProfileLikelihoodCurve()
{
	vdLog2Ratio.clear();
	vdLnLikelihood.clear();
}


bool ProteinReplicate::getProfileLikelihoodCurve(vector<double> & vdLog2RatioRef, vector<double> & vdLnLikelihoodRef )
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

void ProteinReplicate::getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef )
{
	bValidityRef = bValidity;
	dLog2RatioRef = dLog2Ratio;
	dLowerLimitCIRef = dLowerLimitCI;
	dUpperLimitCIRef = dUpperLimitCI;
}

void ProteinReplicate::getEstimates(  double & dLowerLimitCIRef, double & dUpperLimitCIRef )
{
	dLowerLimitCIRef = dLowerLimitCI;
	dUpperLimitCIRef = dUpperLimitCI;
}

string ProteinReplicate::getLocus()
{
	return sLocus;
}


void ProteinReplicate::computeLog2RatioSetting(double & dMLEMinLog2RatioRef, double & dMLEMaxLog2RatioRef, double & dLog2RatioDiscretizationRef)
{

	// normalization
	if( sNormalizationMethod == "Plus" )
	{
		dMLEMinLog2RatioRef = dMLEMinLog2Ratio + dNormalizationValue;
		dMLEMaxLog2RatioRef = dMLEMaxLog2Ratio + dNormalizationValue;		

	}
	else if ( sNormalizationMethod == "Minus" ) 
	{
		dMLEMinLog2RatioRef = dMLEMinLog2Ratio - dNormalizationValue;
		dMLEMaxLog2RatioRef = dMLEMaxLog2Ratio - dNormalizationValue;
	}
	else
	{
		dMLEMinLog2RatioRef = dMLEMinLog2Ratio;
		dMLEMaxLog2RatioRef = dMLEMaxLog2Ratio;
	}	

	// reverse the ratio, if needed
	double dTemp;
	if( bReverseRatio )
	{	
		dTemp = dMLEMinLog2RatioRef;
		dMLEMinLog2RatioRef = -dMLEMaxLog2RatioRef ;
		dMLEMaxLog2RatioRef = -dTemp;		
	}

	dLog2RatioDiscretizationRef = dLog2RatioDiscretization;	
}

/*
double ProteinReplicate::getLog2RatioDiscretization()
{
	return dLog2RatioDiscretization;
}
*/
string ProteinReplicate::getName()
{
	return sName;
}

bool ProteinReplicate::isValid()
{
	return bValidity;
}

int ProteinReplicate::getPeptideNumberWithinCI( double dLowerLimitCI, double dUpperLimitCI )
{
	return pProteinInfo->getValidPeptides(dLowerLimitCI, dUpperLimitCI);
}

int ProteinReplicate::getQuantifiedPeptides()
{
	return iQuantifiedPeptides;
}
	
