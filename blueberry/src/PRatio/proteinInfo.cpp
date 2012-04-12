
#include "proteinInfo.h"


ProteinInfo::ProteinInfo()
{
	bValidity = true;
	sLocus = "";
	sDescription = "";
	dLog2Ratio = 0;
	dLowerLimitCI = ProRataConfig::getMLEMinLog2Ratio();
	dUpperLimitCI = ProRataConfig::getMLEMaxLog2Ratio();
	iQuantifiedPeptides = 0;
	iIdentifiedPeptides = 0;
	iValidPeptides = 0;

	dTotalPeakHeight = 0;
	dTotalPeakArea = 0;
	iMS2SpectralCounts = 0;
}

ProteinInfo::~ProteinInfo()
{
	// destructor
	// vpPeptideInfo is deleted by ProteomeInfo 
}

void ProteinInfo::computeLabelFree()
{
	dTotalPeakHeight = 0;
	dTotalPeakArea = 0;
	iMS2SpectralCounts = 0;
	for( unsigned int i = 0; i < vpPeptideInfo.size(); ++i )
	{
		if( !vpPeptideInfo[i]->getValidity() )
			continue;
		dTotalPeakHeight += vpPeptideInfo[i]->getPeakHeight();
		dTotalPeakArea += vpPeptideInfo[i]->getPeakArea();
		iMS2SpectralCounts += vpPeptideInfo[i]->getMS2Count();
		++iQuantifiedPeptides;
	}
}




bool ProteinInfo::setProteinRatio( ProteinRatio * pProteinRatioInput )
{
	vector< double > vdPeptideLog2Ratio;
	vector< double > vdPeptideLog2SN;

	bValidity = true;

	for( unsigned int i = 0; i < vpPeptideInfo.size(); ++i )
	{
		if( !vpPeptideInfo[i]->getValidity() )
			continue;
/*
		if( vpPeptideInfo[i]->getPCALog2SNR() < ProRataConfig::getMinLog2SNR() )
			continue;

		if( ( vpPeptideInfo[i]->getLocus().size() > 1 ) 
				&& ProRataConfig::getRemoveAmbiguousPeptides() )
			continue;
*/

		vdPeptideLog2Ratio.push_back( vpPeptideInfo[i]->getPCALog2Ratio() );
		vdPeptideLog2SN.push_back( vpPeptideInfo[i]->getPCALog2SNR() );	
	}

	/*
	for( int x = 0; x < vdPeptideLog2Ratio.size(); ++x )
	{
		cout << " peptide ratio = " << vdPeptideLog2Ratio[x] << " peptide SN = " <<  vdPeptideLog2SN[x] << endl;
	}
	*/

	iQuantifiedPeptides = vdPeptideLog2Ratio.size();
	iIdentifiedPeptides = vpPeptideInfo.size();

	if( vdPeptideLog2Ratio.size() < ProRataConfig::getMinPeptideNumber() )
	{
		bValidity = false;
	}
	else
	{

		pProteinRatioInput->process( vdPeptideLog2Ratio, vdPeptideLog2SN );	
		dLog2Ratio = pProteinRatioInput->getProteinLog2Ratio();

		dLowerLimitCI = pProteinRatioInput->getLowerLimitCI();
		dUpperLimitCI = pProteinRatioInput->getUpperLimitCI();

		iValidPeptides = 0;
		for(int j  = 0; j < vdPeptideLog2Ratio.size(); ++j)
		{
			if( ( dLowerLimitCI <= vdPeptideLog2Ratio[j] ) &&  ( vdPeptideLog2Ratio[j] <= dUpperLimitCI ) )
			{
				iValidPeptides = iValidPeptides + 1;
			}	
		}

		if( ( dUpperLimitCI - dLowerLimitCI ) > ProRataConfig::getMaxCIwidth() || iValidPeptides < ProRataConfig::getMinPeptideNumber() )
		{
			bValidity = false;
		}
	}
	return true;
}

void ProteinInfo::addPeptideInfo( PeptideInfo * pPeptideInfoInput )
{
	vpPeptideInfo.push_back( pPeptideInfoInput );
}

void ProteinInfo::setLocus( string sLocusInput )
{
	sLocus = sLocusInput;
}
	
void ProteinInfo::setDescription( string sDescriptionInput)
{
	sDescription = sDescriptionInput;
}

void ProteinInfo::setLog2Ratio( double dLog2RatioInput)
{
	dLog2Ratio = dLog2RatioInput;
}

void ProteinInfo::setLowerLimitCI( double dLowerLimitCIInput)
{
	dLowerLimitCI = dLowerLimitCIInput;
}

void ProteinInfo::setUpperLimitCI( double dUpperLimitCIInput )
{
	dUpperLimitCI = dUpperLimitCIInput;
}

void ProteinInfo::setQuantifiedPeptides( int iQuantifiedPeptidesInput)
{
	iQuantifiedPeptides = iQuantifiedPeptidesInput;
}

void ProteinInfo::setIdentifiedPeptides( int iIdentifiedPeptidesInput)
{
	iIdentifiedPeptides = iIdentifiedPeptidesInput;
}

void ProteinInfo::setValidPeptides( int iValidPeptidesInput)
{
	iValidPeptides = iValidPeptidesInput;
}

vector< PeptideInfo * > & ProteinInfo::getPeptideInfo()
{
	return vpPeptideInfo;
}

vector< string > ProteinInfo::getPeptideSequence()
{
	vector< string > vsPeptideSequence;
	for( unsigned int i = 0; i < vpPeptideInfo.size(); ++i )
	{
		vsPeptideSequence.push_back( vpPeptideInfo[i]->getSequence() );
	}
	return vsPeptideSequence;

}

string ProteinInfo::getLocus()
{
	return sLocus;
}

string ProteinInfo::getDescription()
{
	return sDescription;
}

double ProteinInfo::getLog2Ratio()
{
	return dLog2Ratio;
}
double ProteinInfo::getLowerLimitCI()
{
	return dLowerLimitCI;
}

double ProteinInfo::getUpperLimitCI()
{
	return dUpperLimitCI;
}

int ProteinInfo::getValidPeptides()
{
	return iValidPeptides;
}

int ProteinInfo::getValidPeptides(double dLowerLimitCI, double dUpperLimitCI)
{
	double dPeptideLog2Ratio;
	int iValidPeptideNumber = 0;
	for( int i = 0; i < vpPeptideInfo.size(); ++i )
	{
		if( vpPeptideInfo[i]->getValidity() )
		{
			dPeptideLog2Ratio = vpPeptideInfo[i]->getPCALog2Ratio();
			if( ( dLowerLimitCI <= dPeptideLog2Ratio ) 
					&&  ( dPeptideLog2Ratio <= dUpperLimitCI ) )
			{
				iValidPeptideNumber = iValidPeptideNumber + 1;
			}
		}
	}
	return iValidPeptideNumber;
}

int ProteinInfo::getQuantifiedPeptides()
{
	return iQuantifiedPeptides;
}

int ProteinInfo::getIdentifiedPeptides()
{
	return iIdentifiedPeptides;
}

bool ProteinInfo::getValidity()
{
	return bValidity;
}

bool LessProteinInfo::operator() ( ProteinInfo * pProtein1, ProteinInfo * pProtein2 ) const
{
	if( sKey == "locus" )
	{
		if( pProtein1->getLocus() < pProtein2->getLocus() )
			return true;
		else
			return false;
	}
	else if ( sKey == "log2Ratio" )
	{
		if( pProtein1->getLog2Ratio() < pProtein2->getLog2Ratio() )
			return true;
		else
			return false;
	}
	else if ( sKey == "widthCI" )
	{
		if( ( pProtein1->getUpperLimitCI() - pProtein1->getLowerLimitCI() ) <
				 ( pProtein2->getUpperLimitCI() - pProtein2->getLowerLimitCI() ) )
			return true;
		else
			return false;
	}
	else if ( sKey == "description" )
	{
		if( pProtein1->getDescription() < pProtein2->getDescription() )
			return true;
		else
			return false;
	}
	else if ( sKey == "peptidesQuantified" )
	{
		if( pProtein1->getQuantifiedPeptides() < pProtein2->getQuantifiedPeptides() )
			return true;
		else
			return false;
	}
	else
	{
		if( pProtein1->getLocus() < pProtein2->getLocus() )
			return true;
		else
			return false;
	}
}


