
#include <iostream>
#include "proRataParameters.h"
#include "proteomeInfo.h"
#include <time.h>
using namespace std;

int main( int argc, char * argv[] )
{

	// creat a ProRataParameters object
	ProRataParameters params;
	if( !params.setArguments( argc, argv ) )
	{
		cout << "ERROR: problem in reading the options and configurations!" << endl;
		return 0;
	}

	ProteomeInfo mainProteomeInfo;

	if( !mainProteomeInfo.processPeptidesXIC() )
		cout << "cannot process peptide XIC " << endl;

	if( !mainProteomeInfo.processProteins() )
		cout << "cannot process protein " << endl;

	mainProteomeInfo.writeFileQPR();
	mainProteomeInfo.writeFileTAB();

	cout << "Quantification completed!" << endl;

	return 0;
}
