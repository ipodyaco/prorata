#ifndef CHROMATOGRM_H
#define CHROMATOGRM_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "isotopologue.h"
#include "tinyxml.h"

using namespace std;
class Isotopologue;
class MSdata;
 

/*
 * a class containning the info about an MS2 scan
 * this class corresponds to the <MS2SCORING> element inside <IDENTIFICATION>
 * in the icXML
 */
class MS2Scoring
{
	public:
		MS2Scoring();
		
		// database-searching primary score
		float fScore;

		// MS/MS scan number and retention time 
		unsigned long int iMSMSscan;
		float fRetentionTime;

		// the filename of the identification results
		string sIDfilename;
};

/* 
 * a class containning the info about a protein
 * this class corresponds to the <PROTEIN> element inside <IDENTIFICATION>
 * in the icXML
 */
class Protein
{
	public:
		Protein();
		// the Locus and Description of a protein
		string sLocus;
		string sDescription;
};

/*
 * a class about the identifications of a chromatographic peak
 * one instance of this class will generate a selected ion chromatogram, i.e.
 * a <CHROMATOGRAM> element in icXML
 */
class Identification
{
	public:
		Identification();
		// the indices for the first and last MS2scoring instances in the vMS2Scoring
		// these are not the MS2 scan numbers
		int iFirstMS2;
		int iLastMS2;

		// peptide sequence and charge state
		string sSequence;
		int iChargeState;

		// the adjacent MS2 scans matching this peptide
		// MS2 scans are called adjacent, if they are from the same mzXML/mzData file and
		// are close in retention time, less than ProRataConfig::getMinutesBetweenMS2()
		vector< MS2Scoring > vMS2Scoring;

		// the proteins that this peptides could originate form
		vector< Protein > vProtein;
};

class MZwindows
{
	public:
		vector< float > vfUpperMZ;
		vector< float > vfLowerMZ;
};

class SIC
{
	public:
		string sName;
		MZwindows mzWindows;
		vector< double > vdIntensity;

};

class Chromatogram
{
	public:
		Chromatogram();
		~Chromatogram();

		bool writeChroFile();

		/*
		 * accessors
		 */
		int getIdentifier();
		const Identification & getID();
		float getMaximumScore();
		void getMS2Time( vector< float > & vfMS2Time );
		vector< unsigned long int > getMS2ScanNumber();
		string getSequence();
		int getChargeState();
		bool getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription );
		int getScanCount();
		vector<unsigned long int> getScanVector();
		const vector<float> & getTimeVector();
		const vector<SIC> & getAllSIC();
		bool getIntensityVector( string sName, vector<double> & vdIntensityOutput );
		bool getMZwindows( string sName,  vector< float > & vfLowerMZ, vector< float > & vfUpperMZ  );
		vector<string> getAllSICname();
		vector<string> getAllIDfilename();

		bool isValid();

		unsigned long int getFullScan4Time( float fTime );
		float getFullScanTime4Time( float fTime );
		string getMSfilename();
		
		/*
		 * mutators
		 */
		void setID( const Identification & idInput );
		void setIdentifier( int iIdentifierInput );
		void setMSfilename( string sMSfilenameInput );
		void setScan( const vector< unsigned long int > & viScanInput );
		void setTime( const vector< float > & vfTimeInput );
		void setvSIC( const vector< SIC > & vSICInput );
		void setValidity( bool bValidityInput );

	private:

		Identification myID;

		string sMSfilename;
		int iIdentifier;

		int iScanCount;

		vector< unsigned long int > viScan;

		vector< float > vfTime;	

		vector< SIC > vSIC;
		
		bool bValidity;


};



#endif //CHROMATOGRM_H

