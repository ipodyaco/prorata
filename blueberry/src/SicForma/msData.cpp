#include "msData.h"

MSdata::MSdata()
{
	// constructor
}

MSdata::~MSdata()
{
	// destructor
}

bool MSdata::setFilename( string sFilename)
{
	sMSfilename = sFilename;
	string sCompleteFilenanme = ProRataConfig::getWorkingDirectory() + sFilename;
	if( !myMZreader.setFilename( sCompleteFilenanme ) )
	{
		cout << "ERROR: cannot open the mzXML/mzData file: " << sFilename << endl;
		return false;
	}
	
	unsigned long int iScan = 0;
	int iMSlevel = 0;
	double dPrecurorMZ = 0;
	int iPeaksCount = 0;
	double dRetentionTime = 0;

	vector<int> viAllScanNumbers;
	myMZreader.getAllScanNumbers(viAllScanNumbers);

	for(unsigned int i = 0; i < viAllScanNumbers.size(); ++i )
	{
		iScan = viAllScanNumbers[i];
		if( myMZreader.getHeaderInfo( iScan, &iMSlevel, &dPrecurorMZ, &iPeaksCount, &dRetentionTime ) )
		{
			mTime4Scans[ iScan ] = (float)dRetentionTime;
			if( iMSlevel == 1 )
			{
				mTime4FullScans[ iScan ] = (float)dRetentionTime;
			}
		}
		else
		{
			cout << "WARNING: invalid scan = " << iScan << endl;
		}
	}
	return true;
	
}

string MSdata::getBaseFilename()
{
	string::size_type i = sMSfilename.rfind( ".", ( sMSfilename.length() - 1 ) );  
	if( i != string::npos )
		return sMSfilename.substr( 0, i );
	else
		return "";

}

bool MSdata::getTime4Scan( unsigned long int iScan, float & fTime )
{
	iterTime = mTime4Scans.find( iScan );
	if( iterTime != mTime4Scans.end() )
	{
		fTime = iterTime->second;
		return true;
	}
	else
	{
		// find the two scans immediately before and after iScan
		// interpolate the retention time of the two scans in proportation to the scan numbers
		long int iDistanceCurrent = 0;
		float fTimeBefore = 0;
		float fTimeAfter = 0;
		long int iDistanceBefore = -1000000;
		long int iDistanceAfter = 1000000;
		for( iterTime = mTime4Scans.begin(); iterTime != mTime4Scans.end(); ++iterTime )
		{
			iDistanceCurrent = iterTime->first - iScan;
			if( iDistanceCurrent > iDistanceBefore && iDistanceCurrent < 0 )
			{
				iDistanceBefore = iDistanceCurrent;
				fTimeBefore = iterTime->second;
			}
			if( iDistanceCurrent < iDistanceAfter && iDistanceCurrent > 0 )
			{
				iDistanceAfter = iDistanceCurrent;
				fTimeAfter = iterTime->second;
			}
		}
		float fTimeInterpolationOffset = (fTimeAfter - fTimeBefore ) * (abs(iDistanceBefore)/(abs(iDistanceBefore)+iDistanceAfter));
		fTime = fTimeBefore + fTimeInterpolationOffset;
		return true;
	}
}

void MSdata::getScanVectorTimeVector( float fStartTime, float fEndTime, 
		vector< unsigned long int > & viScanVector, vector< float >  & vfTimeVector )
{
	float fTime;
	for( iterTime = mTime4FullScans.begin(); iterTime != mTime4FullScans.end(); ++iterTime )
	{
		fTime = iterTime->second;
		if( ( fTime >= fStartTime )&&( fTime <= fEndTime ) )
		{
			vfTimeVector.push_back( fTime );
			viScanVector.push_back( iterTime->first );
		}
	}

}

bool MSdata::getIntensityVectors( const vector< unsigned long int > & viScanVector, vector< SIC > & vSIC )
{
	
	vector<float> vfMass;
	vector<float> vfIntensity;
	
	unsigned long int iScan = 0;
	map< unsigned long int, int, less< unsigned long int > >::iterator iterPeakCounts;
	
	unsigned int i = 0;
	unsigned int j = 0;

	double dIntensity;

	bool bNoError = true;
	
	for( i = 0; i < viScanVector.size(); ++i )
	{
		iScan = viScanVector[i];
		if( myMZreader.getPeaksBuffered( iScan, vfMass, vfIntensity ) )
		{
			for( j = 0; j < vSIC.size(); ++j )
			{
				dIntensity = computeIntensity( vSIC[j].mzWindows, vfMass, vfIntensity );
				vSIC[j].vdIntensity.push_back( dIntensity );
			}
		}
		else
		{
			bNoError = false;
			cout << "WARNING: cannot find the full scan = " << iScan << endl;
			for( j = 0; i <  vSIC.size(); ++i )
				vSIC[j].vdIntensity.push_back( 0.0 );

		}
	}
	return bNoError;
}

double MSdata::computeIntensity( const MZwindows & mzWindows, const vector<float> & vfMass, vector<float> & vfIntensity)
{
	float fLowerMZ = 0;
	float fUpperMZ = 0;
	float fPeakMass = 0;
	unsigned int i = 0;
	int iMZWinIndex = 0;
	int iMZwindowNumber = mzWindows.vfLowerMZ.size();
	

	double dSum = 0.0;
	for( i = 0; i < vfMass.size(); ++i )
	{
		fPeakMass = vfMass[i];
		for( iMZWinIndex = 0; iMZWinIndex < iMZwindowNumber; ++iMZWinIndex )
		{
			fUpperMZ = mzWindows.vfUpperMZ[iMZWinIndex];
			fLowerMZ = mzWindows.vfLowerMZ[iMZWinIndex];
			if( ( fPeakMass > fLowerMZ ) && ( fPeakMass < fUpperMZ ) )
			{
				dSum += (double)vfIntensity[i];
				break;
			}
		}
	}
	return dSum; 
}



