
#include "peptideInfo.h"

PeptideInfo::PeptideInfo()
{
	sFilename = "";
	iIdentifier = 0;
	sSequence = "";
	iChargeState = 0;
	fMaximumScore = 0;
	bValidity = false;
	dPCALog2Ratio = 0;
	dPCALog2SNR = 0;
}

PeptideInfo::~PeptideInfo()
{
	// destructor
}

bool PeptideInfo::setPeptideRatio( PeptideRatio * pPeptideRatio )
{
	// Creat a TinyXML document
	TiXmlDocument txdXICFile;

	// Try loading the file.
	string sCompleteFilename = ProRataConfig::getXICxmlDirectory() + sFilename;
	if ( ! ( txdXICFile.LoadFile( sCompleteFilename.c_str() ) ) )
	{
		return false;
		cout << "ERROR! Loading xic.xml file" << sFilename << endl;
	}

	TiXmlElement * pElementChro = txdXICFile.FirstChildElement( "EXTRACTED_ION_CHROMATOGRAMS" );
	int iCurrentIdentifier;
	for( pElementChro = pElementChro->FirstChildElement( "CHROMATOGRAM" ); 
			pElementChro; 
			pElementChro = pElementChro->NextSiblingElement( "CHROMATOGRAM" ) )
	{
		istringstream issStream( pElementChro->Attribute( "identifier" ) );
		issStream >> iCurrentIdentifier;
		if( iCurrentIdentifier == iIdentifier )
		{
			if( !pPeptideRatio->process( pElementChro ) )
			{
				return false;
			}
			else
			{
				setValues( pPeptideRatio );
				return true;
			}
		}
	}	

	cout << "ERROR: cannot find the chromatogram with identifier " << iIdentifier << endl;
	return false;
}

void PeptideInfo::setValues( PeptideRatio * pPeptideRatio )
{
	iIdentifier = pPeptideRatio->getIdentifier();
	sSequence = pPeptideRatio->getSequence();
	iChargeState = pPeptideRatio->getChargeState();
	fMaximumScore = pPeptideRatio->getMaximumScore();
	bValidity = pPeptideRatio->getValidity();
	dPCALog2Ratio = pPeptideRatio->getLog2PCARatio();
	dPCALog2SNR = pPeptideRatio->getLog2PCASN();
	pPeptideRatio->getLocusDescription( vsLocus, vsDescription );
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


