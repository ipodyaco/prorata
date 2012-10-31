
#include <iostream>
#include <string>
#include "projectInfo.h"
#include "proRataConfig.h"
#include <time.h>
using namespace std;

int main( int argc, char * argv[] )
{
	int i = 0;
	string sProRataCombineFilename = "ProRataCombine.xml";
	string sTABOutputFilename = "ProRata_Combine_Measurements.txt";

	string sArgument = "";
	while( i < argc )
	{
		sArgument = argv[i];
		if( sArgument == "-c" )
		{
			i = i + 1;
			sProRataCombineFilename = argv[i];
		}

		if( sArgument == "-o" )
		{
			i = i + 1;
			sTABOutputFilename = argv[i];
		}
		i = i + 1;		
	}

	
	ProjectInfo mainProjectInfo;

	if( !mainProjectInfo.process(sProRataCombineFilename) )
	{
		cout << "cannot process " << sProRataCombineFilename << endl;
	}

	mainProjectInfo.writeFileTAB( sTABOutputFilename );

	cout << "ProRata Combine Done." << endl;
	

	return 0;
}
