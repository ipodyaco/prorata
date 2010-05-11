
#include "proRataConfig.h"

using namespace std;

int main( void )
{
	ProRataConfig::setFilename( "ProRataConfig_test.xml" );

	cout <<  ProRataConfig::getRemoveAmbiguousPeptides() << endl;

	/*
	// test getResidueAtomicComposition()
	map< string, string > map1;
	ProRataConfig::getResidueAtomicComposition( map1 );
	cout << " ResidueAtomicComposition Tables " << endl;
	cout << map1[ "reference" ] << endl;
	cout << map1[ "treatment" ] << endl;
	*/

	/*
	// test getAtomIsotopicComposition()
	vector< double > mass;
	vector< double > natural;
	vector< double > enriched;
	ProRataConfig::getAtomIsotopicComposition( "S", mass, natural, enriched);
	cout << "mass" << '\t' << "natural" << '\t' << "enriched" << endl;
	for( int i = 0; i < mass.size(); i++ )
	{
		cout << mass[i] << '\t' << natural[i] << '\t' << enriched[i] << endl;
	}
	*/

	return 0;
}
