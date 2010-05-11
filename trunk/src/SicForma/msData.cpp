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
	string sCompleteFilenanme = ProRataConfig::getMZxmlDirectory() + sFilename;
	if( !myMZreader.setFilename( sCompleteFilenanme ) )
	{
		cout << "ERROR: cannot open the mzXML/mzData file: " << sFilename << endl;
		return false;
	}
	
	iAnalysisFirstScan = myMZreader.getAnalysisFirstScan();
	iAnalysisLastScan = myMZreader.getAnalysisLastScan();

	unsigned long int iScan = 0;
	int iMSlevel = 0;
	double dPrecurorMZ = 0;
	int iPeaksCount = 0;
	double dRetentionTime = 0;

	for( iScan = iAnalysisFirstScan; iScan <= iAnalysisLastScan; ++iScan )
	{
		if( myMZreader.getHeaderInfo( iScan, &iMSlevel, &dPrecurorMZ, &iPeaksCount, &dRetentionTime ) )
		{
			mTime4Scans[ iScan ] = (float)dRetentionTime;
			if( iMSlevel == 1 )
			{
				mTime4FullScans[ iScan ] = (float)dRetentionTime;
				mPeaksCount4FullScans[ iScan ] = iPeaksCount;
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
		fTime = -999;
		return false;
	}
}

void MSdata::getScanVectorTimeVector( float fStartTime, float fEndTime, 
		vector< unsigned long int > & viScanVector, vector< float >  & vfTimeVector )
{
	float fTime;
	for( iterTime = mTime4FullScans.begin(); iterTime != mTime4FullScans.end(); ++iterTime )
	{
		fTime = iterTime->second;
		if( ( fTime > fStartTime )&&( fTime < fEndTime ) )
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
	int iPeaksCount = 0;
	map< unsigned long int, int, less< unsigned long int > >::iterator iterPeakCounts;
	
	int iMZwindowNumber = vSIC.size();
	int i = 0;
	int j = 0;

	double dIntensity;

	bool bNoError = true;
	
	for( i = 0; i < viScanVector.size(); ++i )
	{
		iScan = viScanVector[i];
		iterPeakCounts = mPeaksCount4FullScans.find( iScan );	
		if( iterPeakCounts != mPeaksCount4FullScans.end() )
		{
			iPeaksCount = iterPeakCounts->second;
			if( myMZreader.getPeaksBuffered( iScan, iPeaksCount, vfMass, vfIntensity ) )
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
				for( j = 0; i <  vSIC.size(); ++i )
					vSIC[j].vdIntensity.push_back( 0.0 );

			}
		}
		else
		{
			cout << "WARNING: cannot find the full scan = " << iScan << endl;
			bNoError = false;
			for( j = 0; j < vSIC.size(); ++j )
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
	int i = 0;
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



