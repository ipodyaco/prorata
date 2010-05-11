
#include <iostream>
using namespace std;

#include "idData.h"
#include "msData.h"
#include "proRataConfig.h"
#include "sicInfo.h"
#include "Isotopologue.h"

int main( int argc, char * argv[] )
{
	ProRataConfig::setFilename( "ProRataConfig.xml" );
	
	vector< Isotopologue * > vpIsotopologue;
	residueMap mAtomicComposition;
	ProRataConfig::getResidueAtomicComposition( mAtomicComposition );
	residueMap::const_iterator iterResidueMap;
	
	for( iterResidueMap = mAtomicComposition.begin(); iterResidueMap != mAtomicComposition.end(); ++iterResidueMap )
	{
		vpIsotopologue.push_back( new Isotopologue( iterResidueMap->first, iterResidueMap->second ) );
	}

	string sSequence;
	cout << "Enter the query amino acid sequence: ";
	cin >> sSequence;
	int iChargeState;
	cout << "Enter the charge state: ";
	cin >> iChargeState;
	
	int i;
	
	string sAtom = "CHONPSCHONPS";
	vector< int > tempAtomicComposition;
	IsotopeDistribution tempIsotopeDistribution;
	MZwindows mzWin;
	vector< double > vdYion;
	vector< double > vdBion;
	for( int n = 0; n < vpIsotopologue.size(); n++ )
	{
		cout << vpIsotopologue[n]->getName() << endl;
		vpIsotopologue[n]->computeAtomicComposition(sSequence, tempAtomicComposition);
		for( i = 0; i < 6; i++ )
			cout << sAtom[i] << " " << tempAtomicComposition[i] << "\t";
		cout << endl;
		for( i = 6; i < tempAtomicComposition.size(); i++ )
			cout << sAtom[i] << " " << tempAtomicComposition[i] << "\t";
		cout << endl;
		vpIsotopologue[n]->computeIsotopicDistribution( sSequence, tempIsotopeDistribution );
		tempIsotopeDistribution.print();
		vpIsotopologue[n]->computeMZwindows( sSequence , iChargeState, mzWin );
		for( i = 0; i < mzWin.vfLowerMZ.size() ; i++ )
		{
			cout << mzWin.vfLowerMZ[i] << " --- " << mzWin.vfUpperMZ[i] << endl;
		}

		cout << "Num	Yion	Bion" << endl;
		if( !vpIsotopologue[n]->computeProductIonMass( sSequence, vdYion, vdBion ) )
			cout << "failed product ion calculation " << endl;
		cout << vdBion.size() << endl;
		for( i = 0 ; i < vdYion.size(); ++i )
		{
			cout << ( i+1 ) << '\t' << vdYion[i] << '\t' << vdBion[i] << endl;
		}
		
	}


	



	return 0;
}
