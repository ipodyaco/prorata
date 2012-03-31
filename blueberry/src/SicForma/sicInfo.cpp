#include "sicInfo.h"

SICinfo::SICinfo()
{
	vpIsotopologue.clear();
	vsMZfilename.clear();
	sIDfilename = "DTASelect-filter.txt";
	
}

SICinfo::~SICinfo()
{
	// destructor
}

bool SICinfo::setFilename( string sIDfilenameInput )
{
	// set ID filename
	sIDfilename =  sIDfilenameInput;

	// get the name of thee working directory 
	string sMZxmlDirectory = ProRataConfig::getWorkingDirectory();

	// remove the last back slash or forward slash in the sMZxmlDirectory
	int iLength = sMZxmlDirectory.length();
	DirectoryStructure dirStructure( sMZxmlDirectory.substr( 0, (iLength - 1) )  );
	
	// get all the filenames matching the extension name specified in ProRataConfig
	// and save them into vsMZfilename
	dirStructure.setPattern( ProRataConfig::getMSfileType() );
	vsMZfilename.clear();
	dirStructure.getFiles( vsMZfilename );

	int i;
	int j;
	// check if the mzFilenames are unique from each other
	// any filename cannot be a substring of any other filename
	// this is because mzFilename is used to group IDs in DTASelect-filter
	for( i = 0; i< vsMZfilename.size(); ++i )
	{
		for( j = (i+1); j < vsMZfilename.size(); ++j )
		{
			if( vsMZfilename[i].find( vsMZfilename[j] ) != string::npos )
			{
				cout << "ERROR: ambiguous mzXML/mzData filenames: " << vsMZfilename[i] << " contains " <<vsMZfilename[j] << endl;
				return false;
			}
			if( vsMZfilename[j].find( vsMZfilename[i] ) != string::npos )
			{
				cout << "ERROR: ambiguous mzXML/mzData filenames: " << vsMZfilename[j] << " contains " <<vsMZfilename[i] << endl;
				return false;
			}

		}	
	}

	return true;
}

bool SICinfo::process(vector< PeptideInfo * > & vpPeptideInfo)
{
	int i;
	int j;
	for( i = 0; i < vpPeptideInfo.size(); i++)
		delete vpPeptideInfo[i];
	vpPeptideInfo.clear();

	// creat the ID data for the given sIDfilename
	cout << "Reading ID file: " << sIDfilename << endl;
	pIDdata = new IDdata;
	if ( !pIDdata->setFilename( sIDfilename ) ){
		return false;
	}

	// the return vector from pIDdata->getIDvector function, which
	// given all IDs from a given MS file
	vector< Identification * > vpIDvector;

	// creat isotopologue objects
	residueMap mAtomicComposition;
	ProRataConfig::getResidueAtomicComposition( mAtomicComposition );
	residueMap::const_iterator iterResidueMap;
	for( iterResidueMap = mAtomicComposition.begin(); iterResidueMap != mAtomicComposition.end(); ++iterResidueMap )
	{
		vpIsotopologue.push_back( new Isotopologue( iterResidueMap->first, iterResidueMap->second ) );
	}

	// for each MS file, a xic.xml file will be created to save all chromatograms extracted from that MS file
	for( i = 0; i < vsMZfilename.size( ); ++i )
	{
		cout << "Processing MS file: " << vsMZfilename[i] << endl;
		// creat the MSdata object
		pMSdata = new MSdata;
		pMSdata->setFilename( vsMZfilename[i] );
		
		// consolidate the IDs from this MS file and save the consolidate IDs
		// into vpIDvector
		pIDdata->consolidateIDlist( pMSdata );
		pIDdata->getIDvector( pMSdata->getBaseFilename(), vpIDvector );

		cout << "Total number of chromatograms to be extracted: " << vpIDvector.size() << endl;
		// iterator thru all IDs in the vpIDvector to extract chromatogram 
		for( j = 0 ; j < vpIDvector.size(); ++j )
		{
			cout << "Extracting chromatogram #" << j << "\r";
			// creat a chromatogram and set all its variables
			Chromatogram chro;
			extractChromatogram( (*vpIDvector[j]), pMSdata->getFilename(), (j+1), chro );
			// calculate peptide ratio and populate vpPeptideInfo
			PeptideRatio currentPeptideRatio;
			if( currentPeptideRatio.process( chro ) )
			{
				PeptideInfo * pCurrentPeptideInfo = new PeptideInfo;
				pCurrentPeptideInfo->setFilename(  chro.getMSfilename() );
				pCurrentPeptideInfo->setValues( &currentPeptideRatio );
				vpPeptideInfo.push_back(pCurrentPeptideInfo);
			}
		}

		delete pMSdata;
	}
	vpIDvector.clear();

	// free memory
	delete pIDdata;
	for( int k = 0; k < vpIsotopologue.size(); ++k)
		delete vpIsotopologue[k];
	return true;
}


bool SICinfo::extractChromatogram( const Identification & idInput, string sMSfilenameInput, int iIdentifier,  Chromatogram & chroOutput )
{
	// set ID and identifer for this chromatogram
	chroOutput.setID( idInput );
	chroOutput.setMSfilename( sMSfilenameInput );
	chroOutput.setIdentifier( iIdentifier );
	
	// set targeted RT window
	float fStartTime;
	float fEndTime;
	fStartTime = idInput.vMS2Scoring[idInput.iFirstMS2].fRetentionTime - ProRataConfig::getMinutesBeforeMS2();
	// start time cannot be negative
	if( fStartTime < 0 )
		fStartTime = 0;
	fEndTime =  idInput.vMS2Scoring[idInput.iLastMS2].fRetentionTime + ProRataConfig::getMinutesAfterMS2();

	// calculate and set the Scan and Time vector
	vector< unsigned long int > viScan;
	vector< float > vfTime;	
	pMSdata->getScanVectorTimeVector( fStartTime, fEndTime, viScan, vfTime);
	chroOutput.setScan( viScan );
	chroOutput.setTime( vfTime );

	// calculate the m/z windows for SIC
	vector< SIC > vSIC;
	for( int i = 0; i < vpIsotopologue.size(); ++i )
	{
		SIC currentSIC;
		currentSIC.sName = vpIsotopologue[i]->getName();
		if(	!vpIsotopologue[i]->computeMZwindows( idInput.sSequence, idInput.iChargeState, currentSIC.mzWindows) )
			return false;
		vSIC.push_back( currentSIC );
	}

	// calculat the intensity vector for SIC
	pMSdata->getIntensityVectors( viScan, vSIC );

	// save the SICs
	chroOutput.setvSIC( vSIC );

	return true;
}

