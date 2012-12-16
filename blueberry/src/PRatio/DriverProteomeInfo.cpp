
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
		  else if (vsArguments[i] == "--chro") { ProRataConfig::setWriteChro(true); }
		  else if (vsArguments[i] == "-h" || vsArguments[i] == "--help") {
			  cout << "Usage: -w WorkingDirectory -c ConfigurationFile -i IdentificationFile --chro" << endl;
			  cout << "ProRata reads MS1 files from WorkingDirectory and writes output to WorkingDirectory" << endl;
			  cout << "If --chro is provided, chromatograms will be written to chro files" << endl;
			  cout << "-w Default: default WorkingDirectory is the current directory; " << endl;
			  cout << "-c Default: default ConfigurationFile is ProRataConfig.xml in WorkingDirectory" << endl;
			  cout << "-i Default: default IdentificationFile is a single pro2psm.txt file in WorkingDirectory" << endl;
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

	sWorkingDirectory = sWorkingDirectory + ProRataConfig::getSeparator();

	if(sConfigFilename == ""){
		sConfigFilename = sWorkingDirectory + "ProRataConfig.xml";
	}

	if(sIDFilename == ""){
		DirectoryStructure dirStructure( sWorkingDirectory.substr( 0, (sWorkingDirectory.length() - 1) )  );
		vector<string> vsPro2PSMfilename;
		dirStructure.setPattern( "pro2psm.txt" );
		dirStructure.getFiles( vsPro2PSMfilename );
		if( vsPro2PSMfilename.size() == 1 ){
			sIDFilename = sWorkingDirectory + vsPro2PSMfilename[0];
		}
		else{
			cout << "Error: please provide an .pro2psm.txt file with the -i option." << endl;
			return 0;
		}
	}
	else{
		// determine if path is provided in the sIDFilename.
		// If no, use working directory.
		size_t found;
		found = sIDFilename.find(ProRataConfig::getSeparator());
		if(found == string::npos){
			sIDFilename = sWorkingDirectory + sIDFilename;
		}
	}

	string sRunBaseName = "ProRata_Quantification";
	size_t separatorFound = sIDFilename.find_last_of(ProRataConfig::getSeparator());
	size_t extensionFound = sIDFilename.find_last_of(".pro2psm.txt");
	if(separatorFound != string::npos && extensionFound != string::npos ){
		// the length of ".pro2psm.txt" is 12
		if( extensionFound - separatorFound > 13  )
		{
			sRunBaseName = sIDFilename.substr(separatorFound + 1, extensionFound - separatorFound - 12 );
		}
	}
//	cout << "sRunBaseName " << sRunBaseName << endl;

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

	if(!ProRataConfig::getIsLabelFree())
	{
		mainProteomeInfo.writeFileQPR(sRunBaseName);
		mainProteomeInfo.writeFileTAB(sRunBaseName);
	}
	else
	{
		mainProteomeInfo.writeFileLabelFree(sRunBaseName);
	}

	cout << "Quantification completed!" << endl;

	return 0;
}
