#include "chromatogram.h"

MS2Scoring::MS2Scoring()
{
	fScore = 0;
	iMSMSscan = 0;
	fRetentionTime = -999;
	sIDfilename = "";
}

Protein::Protein()
{
	sLocus = "";
	sDescription = "";
}

Identification::Identification()
{
	iFirstMS2 = 0;
	iLastMS2 = 0;
	sSequence = "";
	iChargeState = 0;
	
}

Chromatogram::Chromatogram()
{
	iIdentifier = 0;
	sMSfilename = "";
	iScanCount = 0;
	bValidity = true;
	myID.iFirstMS2 = 0;
	myID.iLastMS2 = 0;
	myID.sSequence = "";
	myID.iChargeState = 0;
}

Chromatogram::~Chromatogram()
{

}

int Chromatogram::getIdentifier()
{
	return iIdentifier;
}

const Identification & Chromatogram::getID()
{
	return myID;
}

bool Chromatogram::getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription )
{
	for( unsigned int i = 0; i < myID.vProtein.size(); ++i )
	{
		vsLocus.push_back( myID.vProtein[i].sLocus );
		vsDescription.push_back( myID.vProtein[i].sDescription );
	}
	
	return true;

}

float Chromatogram::getMaximumScore()
{
	float fMaximumScore = 0;
	for( unsigned int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		if( myID.vMS2Scoring[i].fScore > fMaximumScore )
			fMaximumScore = myID.vMS2Scoring[i].fScore;
	}
	return fMaximumScore;

}

void Chromatogram::getMS2Time( vector< float > & vfMS2Time )
{	
	vfMS2Time.clear();
	for( unsigned int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		vfMS2Time.push_back( myID.vMS2Scoring[i].fRetentionTime );
	}
}

vector< unsigned long int > Chromatogram::getMS2ScanNumber()
{	
	vector< unsigned long int > viMS2ScanNumber;
	for( unsigned int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		viMS2ScanNumber.push_back( myID.vMS2Scoring[i].iMSMSscan );
	}
	return viMS2ScanNumber;
}
string Chromatogram::getSequence()
{
	return myID.sSequence;
}

int Chromatogram::getChargeState()
{
	return myID.iChargeState;
}

int Chromatogram::getScanCount()
{
	return iScanCount;
}

vector<unsigned long int> Chromatogram::getScanVector()
{
	return viScan;
}

const vector<float> & Chromatogram::getTimeVector()
{
	return vfTime;
}

bool Chromatogram::getIntensityVector( string sName, vector<double> & vdIntensityOutput )
{
	for( unsigned int i = 0; i < vSIC.size(); ++i )
	{
		if( vSIC[i].sName == sName )
		{
			vdIntensityOutput = vSIC[i].vdIntensity;
			return true;
		}
	}

	cout << "ERROR: cannot find the chromatogram for the isotopologue " << sName << endl;

	return false;
}

bool Chromatogram::getMZwindows( string sName,  vector< float > & vfLowerMZ, vector< float > & vfUpperMZ  )
{
	vfLowerMZ.clear();
	vfUpperMZ.clear();
	for( unsigned int i = 0; i < vSIC.size(); ++i )
	{
		if( vSIC[i].sName == sName )
		{
			vfLowerMZ = vSIC[i].mzWindows.vfLowerMZ;
			vfUpperMZ = vSIC[i].mzWindows.vfUpperMZ;
			return true;
		}
	}

	cout << "ERROR: cannot find the chromatogram for the isotopologue " << sName << endl;

	return false;
}

vector<string> Chromatogram::getAllSICname()
{
	vector< string > vsName;
	for( unsigned int i = 0; i < vSIC.size(); ++i )
	{
		vsName.push_back( vSIC[i].sName );
	}

	return vsName;

}

vector<string> Chromatogram::getAllIDfilename()
{
	vector< string > vsName;
	for( unsigned int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		vsName.push_back( myID.vMS2Scoring[i].sIDfilename );
	}

	return vsName;
}

const vector<SIC> & Chromatogram::getAllSIC()
{
	return vSIC;
}

bool Chromatogram::isValid()
{
	return bValidity;
}


unsigned long int Chromatogram::getFullScan4Time( float fTime )
{
	int iTimeSize = vfTime.size();
	for( int i = 0; i < iTimeSize; ++i )
	{
		if( vfTime[i] > fTime )
		{
			return viScan[i];
		}
	}
	return viScan.back();
}

float Chromatogram::getFullScanTime4Time( float fTime )
{
	int iTimeSize = vfTime.size();
	for( int i = 0; i < iTimeSize; ++i )
	{
		if( vfTime[i] > fTime )
		{
			return vfTime[i];
		}
	}
	return vfTime.back();
}

string Chromatogram::getMSfilename()
{
	return sMSfilename;
}

void Chromatogram::setID( const Identification & idInput )
{ 
	myID = idInput; 
}

void Chromatogram::setIdentifier( int iIdentifierInput )
{ 
	iIdentifier = iIdentifierInput; 
}

void Chromatogram::setMSfilename( string sMSfilenameInput )
{ 
	sMSfilename = sMSfilenameInput; 
}

void Chromatogram::setScan( const vector< unsigned long int > & viScanInput )
{ 
	viScan = viScanInput;
	iScanCount = viScan.size();
}

void Chromatogram::setTime( const vector< float > & vfTimeInput )
{ 
	vfTime = vfTimeInput;
}

void Chromatogram::setvSIC( const vector< SIC > & vSICInput )
{ 
	vSIC = vSICInput; 
}

void Chromatogram::setValidity( bool bValidityInput )
{
	bValidity = bValidityInput;
}

