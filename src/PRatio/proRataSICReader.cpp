
#include "proRataSICReader.h"

ProRataSICReader::ProRataSICReader( const string &sFile )
{
	sFileName = "";
	/*
	if ( !access( sFile.c_str(), R_OK ) )
	{
		cout << "Error: " << sFile << " is not readable." << endl;
		iCurrentIdentifierPointer = 0;
		return ;
	}
	else
	{
		sFileName = sFile;
	}
	*/

	if ( sFileName != "" )
	{
		processHeader( sFileName.c_str() );
	}

	iCurrentIdentifierPointer = 1;
}

ProRataSICReader::ProRataSICReader()
{
	sFileName = "";
	iCurrentIdentifierPointer = 0;
}

ProRataSICReader::~ProRataSICReader()
{

}

int ProRataSICReader::setFileName( const string &sFile )
{
	/*
	if ( !access( sFile.c_str(), R_OK ) )
	{
		cout << "Error: " << sFile << " is not readable." << endl;
		iCurrentIdentifierPointer = 0;
		return  1;
	}
	else
	{
		sFileName = sFile;
	}
	*/
	sFileName = sFile;

	if ( sFileName != "" )
	{
		processHeader( sFileName.c_str() );
	}
	iCurrentIdentifierPointer = 1;

	return 0;
}

string ProRataSICReader::getNextChromatogram()
{
	return getChro( sFileName.c_str(), iCurrentIdentifierPointer );
	iCurrentIdentifierPointer++;

}

string ProRataSICReader::getChromatogram( int iChroIdentifier ) const
{
	int iFirstId = getFirstId();
	int iLastId = getLastId();

	if ( iChroIdentifier < iFirstId ||
			iChroIdentifier > iLastId )
	{
		cout << "Error: Chromatogram identifier " << iChroIdentifier
		       <<  " is not available." << endl;
		return "";
	}

	if ( sFileName == "" )
	{
		cout << "Error: The SIC file is not valid." << endl;
		return "";

	}

	return getChro( sFileName.c_str(), iChroIdentifier );

}

int ProRataSICReader::getFirstIdentifier() const
{	
	return getFirstId();	
}

int ProRataSICReader::getLastIdentifier() const
{	
	return getLastId();	
}

string ProRataSICReader::getMassSpecFile() const
{	
	return getMSFile();	
}

string ProRataSICReader::getProgramName() const
{	
	return getProgram();	
}

string ProRataSICReader::getVersionNumber() const
{	
	return getVersion();	
}
