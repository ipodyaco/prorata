
#include "proRataParameters.h"

ProRataParameters::ProRataParameters()
{
	// default working directory, whihc is the current directory
	sWorkingDirectory = ".";

	// default file name for DTASelect-Filter
	sIDFile = "DTASelect-filter.txt";

	// default file name for ProRataConfig
	sConfigFile = "ProRataConfig.xml";
}

ProRataParameters::~ProRataParameters()
{
}

bool ProRataParameters::setArguments(  int iArgCount,
		char ** cdpArgValues )
{
	int iOption; 

	// Please not print anything if you find unknown parameters.
	opterr = 0;

	// getopt is from the C function, getopt.c
	// getopt.c is a platform-independent implementation of Linux getopt function 
	while( ( iOption = getopt( iArgCount, cdpArgValues, "w:i:") ) != -1 )
	{
		switch (iOption)
		{
			case 'w':
				sWorkingDirectory = optarg;
				break;  
			case 'i':
				sIDFile = optarg;
				break;
			case '?':
				cout << "WARNING: Ignoring unknown input at"
					<< " the command line." << endl;
				break;
		}
	}

	// if the working directory is not given
	// set it to the default, which is assigned in the constructor
	if( !setWorkingDirectory( sWorkingDirectory ) )
		return false;

	// if the ID file name is not given
	// set it to the default, which is assigned in the constructor
	if( !setIDFile( sIDFile ) )
		return false;
		
	// set the configuration file name to the default
	// and read all configurations into memory
	if( !setConfigFile( sConfigFile ) )
		return false;

	return true;
}

bool ProRataParameters::setWorkingDirectory( string sDirectory )
{
	// Save the directory.
	// and append the platform-specific path separators to be ready for
	// appending filename
#ifdef _WIN32
	sWorkingDirectory = sDirectory + "\\";
#else
	sWorkingDirectory = sDirectory + "/";
#endif		
	ProRataConfig::setWorkingDirectory( sWorkingDirectory );
	return true;

	// comment out write permission of the working directory, because in windows
	// the write permission of directories is sometimes not changable.
	/*
	// Check for write permissions of the directory.
	// the access function is changed to  the _access function in Windows
	// with the preprocessors
	if( access( sDirectory.c_str(), W_OK ) != 0 ) {
		cout << "ERROR: Not accessible working directory: " << sDirectory << endl;
		return false;
	}	
	else
	{
		// Save the directory.
		// and append the platform-specific path separators to be ready for
		// appending filename
#ifdef _WIN32
		sWorkingDirectory = sDirectory + "\\";
#else
		sWorkingDirectory = sDirectory + "/";
#endif		
		ProRataConfig::setWorkingDirectory( sWorkingDirectory );
		return true;
	}
	*/
}

bool ProRataParameters::setIDFile( string sFile )
{
	sFile = sWorkingDirectory + sFile;
	
	// Save the ID file.
	sIDFile = sFile;
	return true;
	
}

bool ProRataParameters::setConfigFile( string sFile )
{
	sFile = sWorkingDirectory + sFile;
	// Check for read permission of the file.
	if( access( sFile.c_str(), R_OK ) != 0 ) {
		cout << "ERROR: Not accessible ProRata configuration file: " << sFile << endl;
		return false;
	}	
	
	// Save the ID file.
	sConfigFile = sFile;
	return ProRataConfig::setFilename( sConfigFile );
	
}



