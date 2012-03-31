
#include <time.h>
#include <iostream>
#include "proteomeInfo.h"


using namespace std;

int main( int argc, char * argv[] )
{
	
	  // Grab command line arguments
	  vector<string> vsArguments;
	  while(argc--) vsArguments.push_back(*argv++);

	  // initial values
	  string sWorkingDirectory = "";	  
	  string sConfigFilename = "";
	  string sIDFilename = "";
	  
	  // Parse the arguments
	  unsigned int i  = 0;
	  for(i = 1; i < vsArguments.size(); i++) {
		  if(vsArguments[i] == "-w") { sWorkingDirectory = vsArguments[++i]; }
		  else if (vsArguments[i] == "-c") { sConfigFilename = vsArguments[++i]; }
		  else if (vsArguments[i] == "-i") { sIDFilename = vsArguments[++i]; }
		  else if (vsArguments[i] == "-h" || vsArguments[i] == "--help") {
			  cout << "Usage: -w WorkingDirectory -c ConfigurationFile -i IdentificationFile" << endl;
			  cout << "Default: WorkingDirectory is the current directory; " << endl;
			  cout << "ConfigurationFile is ProRataConfig.xml in the working directory" << endl;
			  cout << "IdentificationFile is DTASelect-filter.txt in the working directory" << endl;
			  exit(0);
	    }
	    else { 
		    // unknown arguments
		 cerr << "Unknown option " << vsArguments[i] << endl << endl;
		 exit(1);
	    }
	  }

	if(sWorkingDirectory == ""){
		sWorkingDirectory = ".";
	}
#ifdef _WIN32
	sWorkingDirectory = sWorkingDirectory + "\\";
#else
	sWorkingDirectory = sWorkingDirectory + "/";
#endif

	if(sConfigFilename == ""){
		sConfigFilename = sWorkingDirectory + "ProRataConfig.xml";
	}

	if(sIDFilename == ""){
		sIDFilename = sWorkingDirectory + "DTASelect-filter.txt";
	}

	// Load configuration file.
	cout << "Reading config file: " << sConfigFilename << endl;
	if(!ProRataConfig::setFilename( sConfigFilename )) {
		cerr << "Could not load config file " << sConfigFilename << endl << endl;
		cout << "Usage: -w WorkingDirectory -c ConfigurationFile -i IdentificationFile" << endl;
		exit(2);
	}	

	ProRataConfig::setWorkingDirectory( sWorkingDirectory );

	ProteomeInfo mainProteomeInfo;

	if( !mainProteomeInfo.processPeptidesXIC( sIDFilename ) )
		cout << "Error: cannot process peptide XIC " << endl;

	if( !mainProteomeInfo.processProteins() )
		cout << "Error: cannot process protein " << endl;

	mainProteomeInfo.writeFileQPR();
	mainProteomeInfo.writeFileTAB();

	cout << "Quantification completed!" << endl;

	return 0;
}
