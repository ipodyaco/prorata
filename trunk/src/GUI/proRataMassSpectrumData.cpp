
#include "proRataMassSpectrumData.h"
#include <QMessageBox>
#include <QtGlobal>

#include <iostream>
#include <QVector>
using namespace std;

ProRataMassSpectrumData::ProRataMassSpectrumData()
{

	// get the name of thee working directory 
	string sMZxmlDirectory = ProRataConfig::getMZxmlDirectory();

	// remove the last back slash or forward slash in the sWorkingDirectory
	int iLength = sMZxmlDirectory.length();
	DirectoryStructure dirStructure( sMZxmlDirectory.substr( 0, (iLength - 1) )  );
	
	// get all the filenames matching the extension name specified in ProRataConfig
	// and save them into vsMZfilename
	dirStructure.setPattern( ProRataConfig::getMSfileType() );
	vector<string> vsMZfilename;
	dirStructure.getFiles( vsMZfilename );

	for( unsigned int iC = 0; iC < vsMZfilename.size(); iC++ )
	{
		mzReader *mzrTemp = new mzReader;
		mzrTemp->setFilename( sMZxmlDirectory + vsMZfilename.at(iC) );
		mapReaderFileMappings[vsMZfilename.at(iC)] = mzrTemp;
	}

}

ProRataMassSpectrumData::~ProRataMassSpectrumData()
{
	map<string, mzReader*>::iterator itr;
	itr = mapReaderFileMappings.begin();

	for( ; itr != mapReaderFileMappings.end(); itr++ )
	{
		delete (*itr).second;
	}
}

bool ProRataMassSpectrumData::getScan( string sFilename, unsigned long int iScan, 
		vector< double > & vdMZ, vector< double > & vdRelativeIntensity )
{
	vdMZ.clear();
	vdRelativeIntensity.clear();
	
	mzReader * mzRdr = mapReaderFileMappings[sFilename];

	if ( !mzRdr )
	{
		return false;
	}

	int iMSLevel, iPeaksCount;
	double dPrecursorMZ, dRententionTime;

	if ( !( mzRdr->getHeaderInfo( iScan, &iMSLevel, 
				&dPrecursorMZ, 	&iPeaksCount,
				&dRententionTime ) ))
	{
		return false;
	}

	vector<float> vfMass;
	vector<float> vfInten;

	if ( !( mzRdr->getPeaks( iScan, iPeaksCount, vfMass, vfInten ) ))
	{
		return false;
	}
	
	double dMaxInten = (double)( *max_element( vfInten.begin(), vfInten.end() ) );

	// cast float to double
	for( unsigned int i = 0; i < vfMass.size(); i++ )
	{
		vdMZ.push_back( (double)(vfMass.at(i)) );
	}

	if( dMaxInten > 0.000001 )
	{
		// cast float to double and calculate relative intensity
		for( unsigned int i = 0; i < vfInten.size(); i++ )
		{
			vdRelativeIntensity.push_back( ( (double)(vfInten.at(i)) / dMaxInten ) * 100.0 );
		}
	}
	else
	{
		for( unsigned int i = 0; i < vfInten.size(); i++ )
		{
			vdRelativeIntensity.push_back( 0.0 );
		}
	}

	return true;
}

/*
 * overloaded getScan to provide dPrecursorMZ,
 * the code is mostly cut-and-paste from the previous function
 */

bool ProRataMassSpectrumData::getScan( string sFilename, unsigned long int iScan, 
		vector< double > & vdMZ, vector< double > & vdRelativeIntensity, double & dPrecursorMZ )
{
	vdMZ.clear();
	vdRelativeIntensity.clear();
	
	mzReader * mzRdr = mapReaderFileMappings[sFilename];

	if ( !mzRdr )
	{
		return false;
	}

	int iMSLevel, iPeaksCount;
	double dRententionTime;

	if ( !( mzRdr->getHeaderInfo( iScan, &iMSLevel, 
				&dPrecursorMZ, 	&iPeaksCount,
				&dRententionTime ) ))
	{
		return false;
	}

	vector<float> vfMass;
	vector<float> vfInten;

	if ( !( mzRdr->getPeaks( iScan, iPeaksCount, vfMass, vfInten ) ))
	{
		return false;
	}
	
	double dMaxInten = (double)( *max_element( vfInten.begin(), vfInten.end() ) );

	// cast float to double
	for( unsigned int i = 0; i < vfMass.size(); i++ )
	{
		vdMZ.push_back( (double)(vfMass.at(i)) );
	}
	
	// cast float to double and calculate relative intensity
	for( unsigned int i = 0; i < vfInten.size(); i++ )
	{
		vdRelativeIntensity.push_back( ( (double)(vfInten.at(i)) / dMaxInten ) * 100.0 );
	}

	return true;
}

