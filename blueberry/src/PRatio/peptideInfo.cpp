
#include "peptideInfo.h"

PeptideInfo::PeptideInfo()
{
	sFilename = "";
	iIdentifier = 0;
	sSequence = "";
	iChargeState = 0;
	fMaximumScore = 0;
	bValidity = false;
	iMS2Count = 0;
	fPeakTimeWidth = 0;
	fLeftValleyTime = 0;
	fRightValleyTime = 0;

	dPCALog2Ratio = 0;
	dPCALog2SNR = 0;

	dPeakHeight = 0;
	dPeakArea = 0;
	dPeakSNR = 0;

}

PeptideInfo::~PeptideInfo()
{
	// destructor
}

void PeptideInfo::setValues( PeptideRatio * pPeptideRatio )
{
	// common between PeptideRatio and PeptideLabelFree
	iIdentifier = pPeptideRatio->getIdentifier();
	sSequence = pPeptideRatio->getSequence();
	iChargeState = pPeptideRatio->getChargeState();
	fMaximumScore = pPeptideRatio->getMaximumScore();
	bValidity = pPeptideRatio->getValidity();
	pPeptideRatio->getLocusDescription( vsLocus, vsDescription );
	iMS2Count = pPeptideRatio->getMS2ScanNumber().size();
	fLeftValleyTime = pPeptideRatio->getLeftValleyTime();
	fRightValleyTime = pPeptideRatio->getRightValleyTime();
	fPeakTimeWidth = fRightValleyTime - fLeftValleyTime;
	vfMS2Time = pPeptideRatio->getMS2Time(); 
	vsAllIDfilename = pPeptideRatio->getAllIDfilename();

	// specific for PeptideRatio
	dPCALog2Ratio = pPeptideRatio->getLog2PCARatio();
	dPCALog2SNR = pPeptideRatio->getLog2PCASN();
}

void PeptideInfo::setValues( PeptideLabelFree * pPeptideLabelFree )
{
	// common between PeptideRatio and PeptideLabelFree	
	iIdentifier = pPeptideLabelFree->getIdentifier();
	sSequence = pPeptideLabelFree->getSequence();
	iChargeState = pPeptideLabelFree->getChargeState();
	fMaximumScore = pPeptideLabelFree->getMaximumScore();
	bValidity = pPeptideLabelFree->getValidity();
	pPeptideLabelFree->getLocusDescription( vsLocus, vsDescription );
	iMS2Count = pPeptideLabelFree->getMS2ScanNumber().size();
	fLeftValleyTime = pPeptideLabelFree->getLeftValleyTime();
	fRightValleyTime = pPeptideLabelFree->getRightValleyTime();
	fPeakTimeWidth = fRightValleyTime - fLeftValleyTime;
	vfMS2Time = pPeptideLabelFree->getMS2Time(); 
	vsAllIDfilename = pPeptideLabelFree->getAllIDfilename();

	// specific for PeptideLabelFree	
	dPeakHeight = pPeptideLabelFree->getPeakHeight();
	dPeakArea = pPeptideLabelFree->getPeakArea();
	dPeakSNR = pPeptideLabelFree->getPeakSNR();
}

void PeptideInfo::setFilename( string sFilenameInput )
{
	sFilename = sFilenameInput;
}

void PeptideInfo::setIdentifier( int iIdentifierInput )
{
	iIdentifier = iIdentifierInput;
}

void PeptideInfo::setSequence( string sSequenceInput )
{
	sSequence = sSequenceInput;
}

void PeptideInfo::setChargeState( int iChargeStateInput )
{
	iChargeState = iChargeStateInput;
}

void PeptideInfo::setMaximumScore( float fScoreInput )
{
	fMaximumScore = fScoreInput;
}

void PeptideInfo::setValidity( bool bValidityInput )
{
	bValidity = bValidityInput;
}

void PeptideInfo::setPCALog2Ratio( double dPCALog2RatioInput )
{
	dPCALog2Ratio = dPCALog2RatioInput;

}

void PeptideInfo::setPCALog2SNR( double dPCALog2SNRInput )
{
	dPCALog2SNR = dPCALog2SNRInput;
}

void PeptideInfo::setLocus( vector< string > vsLocusInput )
{
	vsLocus = vsLocusInput;
}

string PeptideInfo::getFilename()
{
	return sFilename;
}

int PeptideInfo::getIdentifier()
{
	return iIdentifier;
}

string PeptideInfo::getSequence()
{
	return sSequence;
}
	
int PeptideInfo::getChargeState()
{
	return iChargeState;
}
	
float PeptideInfo::getMaximumScore()
{
	return fMaximumScore;
}
	
bool PeptideInfo::getValidity()
{
	return bValidity;
}
	
double PeptideInfo::getPCALog2Ratio()
{
	return dPCALog2Ratio;
}
	
double PeptideInfo::getPCALog2SNR()
{
	return dPCALog2SNR;
}
	
const vector< string > & PeptideInfo::getLocus()
{
	return vsLocus;
}

const vector< string > & PeptideInfo::getDescription()
{
	return vsDescription;
}

bool LessPeptideInfo::operator() ( PeptideInfo * pPeptide1, PeptideInfo * pPeptide2 ) const
{
	if( sKey == "sequence" )
	{
		if( pPeptide1->getSequence() < pPeptide2->getSequence() )
			return true;
		else
			return false;
	}
	else if ( sKey == "log2Ratio" )
	{
		if( pPeptide1->getPCALog2Ratio() < pPeptide2->getPCALog2Ratio() )
			return true;
		else
			return false;
	}
	else if ( sKey == "log2SN" )
	{
		if( pPeptide1->getPCALog2SNR() < pPeptide2->getPCALog2SNR() )
			return true;
		else
			return false;
	}
	else if ( sKey == "validity" )
	{
		if( pPeptide1->getValidity() < pPeptide2->getValidity() )
			return true;
		else
			return false;
	}
	else
	{
		cout << "ERROR: The key cannot be recoginzed and peptides are sorted by sequence ! " << endl;
		if( pPeptide1->getSequence() < pPeptide2->getSequence() )
			return true;
		else
			return false;
	}
	
}


