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
	string sMZxmlDirectory = ProRataConfig::getMZxmlDirectory();

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

bool SICinfo::process()
{

	// creat the directory for xic files
	string sXICxmlDirectory =  ProRataConfig::getXICxmlDirectory();
	int iLength = sXICxmlDirectory.length();
	string sTempDir = sXICxmlDirectory.substr(0, (iLength - 1));
	cout << " making XIC directory " << mkdir( sTempDir.c_str() ) << "  " << sTempDir << endl;
/*	
	if(  mkdir( sTempDir.c_str() ) != 0 )
	{
		cout << "ERROR: cannot creat xic directory: " << sTempDir << endl;
	       return false;
	}
*/	
	
	// creat the ID data for the given sIDfilename
	cout << " reading the ID file  " << sIDfilename << endl;
	pIDdata = new IDdata;
	pIDdata->setFilename( sIDfilename );

	// the return vector from pIDdata->getIDvector function, which
	// given all IDs from a given MS file
	vector< Identification * > vpIDvector;
	
	vector< MZwindows > vMZwin;

	// the EXTRACTED_ION_CHROMATOGRAMS element's end tag in the output file
	string sRootElementEndTag = "</EXTRACTED_ION_CHROMATOGRAMS>\n";

	int i;
	int j;
	int n;
	
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
		cout << "reading MS files " << vsMZfilename[i] << endl;
		// creat the MSdata object
		pMSdata = new MSdata;
		pMSdata->setFilename( vsMZfilename[i] );
		
		// consolidate the IDs from this MS file and save the consolidate IDs
		// into vpIDvector
		pIDdata->consolidateIDlist( pMSdata );
		pIDdata->getIDvector( pMSdata->getBaseFilename(), vpIDvector );

		// calculate how many xic files are needed
		int iXICfileCount = 0;
		if( (vpIDvector.size() % CHRO_COUNT) != 0 )
		{
			iXICfileCount = (int)( vpIDvector.size() / CHRO_COUNT ) + 1;
		}
		else
		{
			iXICfileCount = (int)( vpIDvector.size() / CHRO_COUNT );
		}
		
		for( n = 0 ; n < iXICfileCount; ++n )
		{
			int iFirstIdentifier = n * CHRO_COUNT + 1;
			int iLastIdentifier = n * CHRO_COUNT + CHRO_COUNT;
			if( iLastIdentifier > vpIDvector.size() )
				iLastIdentifier = vpIDvector.size();
			
			// set the filename for the output xic.xml file
			ostringstream ossStream;
			ossStream << pMSdata->getBaseFilename() << "." << iFirstIdentifier << "." << iLastIdentifier << ".xic.xml";
			string xicFilename =  ProRataConfig::getXICxmlDirectory() + ossStream.str();
			// write the start tag for the root element into the xic.xml file
			// if the file already exists, its current content will be overwritten
			// otherwise the file will be created
			if( !writeRootElementStartTag( xicFilename, iFirstIdentifier, iLastIdentifier ) )
				return false;

			// re-open the file in the "append" mode, such that all chromatograms will
			// be appended to the end of the file
			FILE * pFile;
			if( ( pFile = fopen( xicFilename.c_str(), "a" ) ) == NULL  ) 
			{
				cout << "ERROR: cannot write file: " << xicFilename << endl;
				return false;
			}

			cout << "Formulating XIC: " << ossStream.str() << endl;
			
			// iterator thru all IDs in the vpIDvector to extract chromatogram 
			for( j = ( iFirstIdentifier - 1 ); j < iLastIdentifier; ++j )
			{

				// creat a chromatogram and set all its variables
				Chromatogram chro;
				extractChromatogram( (*vpIDvector[j]), pMSdata->getFilename(), (j+1), chro );

				// append its content into the file
				if( !chro.writeToXicFile( pFile ) )
					return false;
			//	fputs( "\n", pFile );
				fflush( pFile );
			}

			// write the end tag for EXTRACTED_ION_CHROMATOGRAMS
			fputs( sRootElementEndTag.c_str(), pFile );
			// clean up
			fclose( pFile );
		}
		vpIDvector.clear();
		delete pMSdata;
	}

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

bool SICinfo::writeRootElementStartTag( string sXicFilename, int iFirstIdentifier, int iLastIdentifier )
{
	// open the file in the "write" mode
	FILE * pFile;
	if( ( pFile = fopen( sXicFilename.c_str(), "w" ) ) == NULL  ) 
	{
		cout << "ERROR: cannot write file: " << sXicFilename << endl;
		return false;
	}

	// write the declaration
	string sDeclaration = "<?xml version = \"1.0\" ?>\n";
	fputs( sDeclaration.c_str(), pFile );

	// formulate a string for the start end with all its attributes
	ostringstream ossStream;
	ossStream << "<EXTRACTED_ION_CHROMATOGRAMS MSfile=\"" << pMSdata->getFilename() <<"\" ";
	ossStream << "first_identifier=\"" << iFirstIdentifier << "\" ";
	ossStream << "last_identifier=\"" << iLastIdentifier << "\" ";
	ossStream << "program=\"ProRata\" version=\"" << ProRataConfig::getProRataVersion() << "\"> \n" ;

	// write the start tag into the file
	string sRootElementStartTag = ossStream.str();
	fputs( sRootElementStartTag.c_str(), pFile );
	fclose( pFile );

	if( !ProRataConfig::writeConfigXML( sXicFilename, 1, true, false, false ) )
	{
		cout << "ERROR: cannot write the CONFIG element to the xic file: " << sXicFilename << endl;
		return false;
	}
	
	return true;
}

