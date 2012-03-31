#ifndef ISOTOPOLOGUE_H
#define ISOTOPOLOGUE_H

#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include "chromatogram.h"
#include "msData.h"

using namespace std;

class MZwindows;

class IsotopeDistribution
{
	public:
		IsotopeDistribution();
		IsotopeDistribution( vector< double > vItsMass,  vector< double > vItsComposition );
		~IsotopeDistribution();
		
		vector< double > vMass;
		vector< double > vProb;

		double getMostAbundantMass();
		double getAverageMass();

		// print out the isotoptic distribution
		// this is mainly used for debuging
		void print();


};

class Isotopologue
{
	public:
		Isotopologue();
		Isotopologue( string sName, string sTable );
		~Isotopologue();
		
		// setup all variables from configuration
		bool setupIsotopologue(  const string & sName, const string & sTable );
		
		// compute MZ window for an amino acid sequence
		bool computeMZwindows( string sSequence, int iChargeState, MZwindows & myMZwindows );

		// compute the mass of the most abundant isotopologue
		double computeMostAbundantMass( string sSequence );
		double computeAverageMass( string sSequence );

		// compute the mass of a peptide's Y-ion series
		bool computeProductIonMass( string sSequence, vector< double > & vdYion, vector< double > & vdBion );

		// compute isotoptic distribution for an amino acid sequence
		bool computeIsotopicDistribution( string sSequence, IsotopeDistribution & myIsotopeDistribution );
		
		// compute isotoptic distribution for a given atomic composition, which can be that of a residue's or a amino acid sequence's
		bool computeIsotopicDistribution( vector< int > AtomicComposition, IsotopeDistribution & myIsotopeDistribution );
		
		// compute the atomic composition for an amino acid sequence
		bool computeAtomicComposition( string sSequence, vector< int >  & myAtomicComposition );
		
		string getName()
		{ return sIsotoplogueName; }

	private:
		// functions for IsotopeDistribution's arithmetic
		IsotopeDistribution sum( IsotopeDistribution distribution0, IsotopeDistribution distribution1); 
		IsotopeDistribution multiply( IsotopeDistribution distribution0, int count );

		void formatSequence( string & sSequence );

		// implementation of max and min
		double maximum( double a, double b )
		{ return (a > b) ? a : b; }
		double minimum( double a, double b)
		{ return (a < b) ? a : b; }
		
		// when two peaks have a mass difference less than the MassPrecision
		// they will be merged into one peak with their average mass and sum intensity
		const double MassPrecision; 
		
		// when a peak have a probability less than the ProbabilityCutoff
		// this peak will be ingnored, which makes the total probability space less than 1
		const double ProbabilityCutoff; // 1*10E-9
		
		// the name of atoms
		const string AtomName;
		
		// the number of natural CHONPS and enriched CHONPS
		const int AtomNumber;
		
		// variables for this isotopologue
		string sIsotoplogueName;
		map< string, vector< int > > mResidueAtomicComposition;
		vector< IsotopeDistribution > vAtomIsotopicDistribution;
		map< string, IsotopeDistribution > vResidueIsotopicDistribution;
		
};

#endif //ISOTOPOLOGUE_H
