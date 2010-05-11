
#include "proRataSICReader.h"

#include <iostream>
using namespace std;

int main()
{
	ProRataSICReader prsrFirst;
	prsrFirst.setFileName( "sample.xml" );

	cout << prsrFirst.getMassSpecFile() << endl;

	cout << prsrFirst.getChromatogram( 2 ) << endl;

	return 0;
}
