
#ifndef PRORATAPARAMETERS_H
#define PRORATAPARAMETERS_H

#include "proRataConfig.h"

#ifdef _WIN32
#include "getopt.h"
#include <stdio.h>
#include <stdlib.h>
#include <io.h>
//#define W_OK 02
//#define R_OK 04
#define access( a, b ) _access( a, b )
#else
#include <unistd.h>
#endif

#include <string>

// External identifiers for getopt
extern char *optarg;
extern int  opterr;

using namespace std;

class ProRataParameters
{
	public:
		ProRataParameters();
		virtual ~ProRataParameters();

		// set sWorkingDirectory and sIDFile with argc, argv
		bool virtual setArguments( int , char ** );

		bool virtual setWorkingDirectory( string sDirectory );
		bool virtual setIDFile( string sFile );
		bool virtual setConfigFile( string sFile );

		string getWorkingDirectory() const
		{	return sWorkingDirectory;	}
		string getIDFilename() const
		{	return sIDFile;		}
		string getConfigFilename() const
		{	return sConfigFile;	}

	protected:

		string sWorkingDirectory;
		string sIDFile;
		string sConfigFile;
		
};

#endif //PRORATAPARAMETERS_H
