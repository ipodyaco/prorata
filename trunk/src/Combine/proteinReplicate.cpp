#include "proteinReplicate.h"

ProteinReplicate::ProteinReplicate()
{
	sName = "";
	pProteomeInfo = NULL;

	bRemoveAmbiguousPeptide = true;
	iMinPeptideNumber = 0;
	dMaxCIwidth = 0;
	dMinLog2SNR = 0;
	dMaxLog2SNR = 0;
	dMLEMinLog2Ratio = 0;
	dMLEMaxLog2Ratio = 0;
	dLog2RatioDiscretization = 0;
	dSDSlope = 0;
	dSDIntercept = 0;
	dMeanSlope = 0;
	dMeanIntercept = 0; 
	dSmoothingProbSpace = 0;
	dLnLikelihoodCutoffOffset = 0;
	
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

void ProteinReplicate::setProteomeInfo( ProteomeInfo * pProteomeInfoConfig )
{
	pProteomeInfo = pProteomeInfoConfig;

	bRemoveAmbiguousPeptide	= ProRataConfig::getRemoveAmbiguousPeptides();		
	iMinPeptideNumber	= ProRataConfig::getMinPeptideNumber();			
	dMaxCIwidth		= ProRataConfig::getMaxCIwidth();				
	dMinLog2SNR		= ProRataConfig::getMinLog2SNR();				
	dMaxLog2SNR		= ProRataConfig::getMaxLog2SNR();				
	dMLEMinLog2Ratio	= ProRataConfig::getMLEMinLog2Ratio();			
	dMLEMaxLog2Ratio	= ProRataConfig::getMLEMaxLog2Ratio();			
	dLog2RatioDiscretization	= ProRataConfig::getLog2RatioDiscretization();		
	dSDSlope		= ProRataConfig::getSDSlope();				
	dSDIntercept		= ProRataConfig::getSDIntercept();			
	dMeanSlope		= ProRataConfig::getMeanSlope();				
	dMeanIntercept		= ProRataConfig::getMeanIntercept();			
	dSmoothingProbSpace	= ProRataConfig::getSmoothingProbilitySpace();		
	dLnLikelihoodCutoffOffset	= ProRataConfig::getLnLikelihoodCutoffOffset();		

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

	

	// validity is retrieved from pProteinInfo read from QPR file
	bValidity = pProteinInfo->getValidity();
	iQuantifiedPeptides = pProteinInfo->getQuantifiedPeptides();
	delete pProteinRatio;
	
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

void ProteinReplicate::getLog2RatioSetting(double & dMLEMinLog2RatioRef, double & dMLEMaxLog2RatioRef, double & dLog2RatioDiscretizationRef)
{
	dMLEMinLog2RatioRef = dMLEMinLog2Ratio;
	dMLEMaxLog2RatioRef = dMLEMaxLog2Ratio;
	dLog2RatioDiscretizationRef = dLog2RatioDiscretization;	
}

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
	
