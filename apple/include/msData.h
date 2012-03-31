#ifndef MSDATA_H
#define MSDATA_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include "chromatogram.h"
#include "mzReader.h"
#include "proRataConfig.h"

using namespace std;

class SIC;
class MZwindows;

class MSdata
{
	public:
		MSdata();
		~MSdata();
		
		bool setFilename( string sFilename);
		string getFilename()
		{ return sMSfilename; }

		string getBaseFilename();
		
		// get the retention time for scans
		bool getTime4Scan( unsigned long int iScan, float & fTime );
		
		// get all full scan numbers between two retention time points
		void getScanVectorTimeVector( float fStartTime, float fEndTime, 
				vector< unsigned long int > & viScanVector, vector< float >  & vfTimeVector );

		// get the intensties for a selected ion chromatogram
		bool getIntensityVectors( const vector< unsigned long int > & viScanVector,
				vector< SIC > & vSIC );
		
	private:
		string sMSfilename;

		
		// the mzReader for random access of mzXML/mzData
		mzReader myMZreader;

		// the look-up table of the retention time for full scans
		map< unsigned long int, float, less< unsigned long int > > mTime4FullScans;

		// the look-up table of the retention time for all scans
		map< unsigned long int, float, less< unsigned long int > > mTime4Scans;

		map< unsigned long int, float, less< unsigned long int > >::iterator iterTime;

		// compute the intensity within a set of MZ windows for the given mass spectrum 
		double computeIntensity( const MZwindows & mzWindows, const vector<float> & vfMass, vector<float> & vfIntensity);
		
};

#endif //MSDATA_H
