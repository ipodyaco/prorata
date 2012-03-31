
#include <iostream>
using namespace std;

#include "idData.h"
#include "msData.h"
#include "proRataParameters.h"
#include "sicInfo.h"
#include "isotopologue.h"

int main( int argc, char * argv[] )
{

	// creat a ProRataParameters object
	ProRataParameters params;
	if( !params.setArguments( argc, argv ) )
	{
		cout << "ERROR: problem in reading the options and configurations!" << endl;
		return 0;
	}


	SICinfo sicInfo;
	sicInfo.setFilename( params.getIDFilename() );

	sicInfo.process();

	cout << "Ion chromatogrom extraction completed!" << endl;


	return 0;
}
