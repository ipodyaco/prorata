#ifndef PROTEINCOMBINED_H
#define PROTEINCOMBINED_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <iomanip>
#include <fstream>

#include "proteinDirectComparison.h"
#include "proteinIndirectComparison.h"
#include "proteinReplicate.h"
#include "proteinInfo.h"
#include "proteinRatio.h"
#include "directoryStructure.h"
#include "proRataConfig.h"
#include "tinyxml.h"

using namespace std;

class ProteinCombined
{
	
	public:
		ProteinCombined();
		~ProteinCombined();

		// copy constructor, used to creat new ProteinCombined from templateProteinCombined in ProjectInfo
		ProteinCombined( const ProteinCombined & templateProteinCombined );

		bool runQuantification( string sInputLocus, string sInputDescription );
		bool runQuantification( string sInputLocus );
		void eraseProfileLikelihoodCurve();

		string getLocus();
		string getDescription();
		bool getValidity();
		const vector< ProteinDirectComparison > & getProteinDirectComparisonVector();		
		const vector< ProteinIndirectComparison > & getProteinIndirectComparisonVector();		
		
		// return false, if cannot find the comparison matching the name;
		bool getDirectComparison( string sComparisonName, ProteinDirectComparison & comparison );
		bool getIndirectComparison( string sComparisonName, ProteinIndirectComparison & comparison );
		bool IsThereDirectComparison( string sComparisonName );
		bool IsThereIndirectComparison( string sComparisonName );

		/*
		 *  these functions are supposed to be applicable only to templateProteinCombined in ProjectInfo;
		 *  addDirectComparison:	push back a ProteinDirectComparison to vDirectComparison
		 *  addIndirectComparison:	push back a ProteinIndirectComparison to vIndirectComparison
		 */
		bool addDirectComparison(ProteinDirectComparison tempProteinDirectComparison);
		bool addIndirectComparison(ProteinIndirectComparison tempProteinIndirectComparison);
		
	private:

		void computeValidity();
		
		bool bValidity;
		string sLocus;
		string sDescription;

		vector< ProteinDirectComparison > vDirectComparison;
		vector< ProteinIndirectComparison > vIndirectComparison;

};

class LessProteinCombined
{
	public:
		LessProteinCombined();
		LessProteinCombined( string sKeyInput ) { sKey = sKeyInput; }
		
		void setKey( string sKeyInput ) { sKey = sKeyInput; }
		string getKey() { return sKey; }
		
		/*
		 * the proteinCombined pointers can be sorted by a number of member variables:
		 * the sKey specifies which variables to be used to sort
		 * the accepted values for sKay and the corresponding accessors of the variables are
		 * "locus"			pProtein1->getLocus()
		 * "description"		pProtein1->getDescription() 
		 */
		bool operator() ( ProteinCombined * pProtein1, ProteinCombined * pProtein2 ) const;
		
	private:
		string sKey;
};


#endif //PROTEINCOMBINED_H
