
#ifndef PROJECTINFO_H
#define PROJECTINFO_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <iomanip>
#include <fstream>

#include "proteinCombined.h"
#include "proteinDirectComparison.h"
#include "proteinIndirectComparison.h"
#include "proteinReplicate.h"
#include "proteinInfo.h"
#include "proteomeInfo.h"
#include "proteinRatio.h"
#include "proRataConfig.h"
#include "tinyxml.h"

using namespace std;

class ProjectInfo
{
	public:
		ProjectInfo();
		~ProjectInfo();

		bool process( string sProRataCombineXMLfilename );
		
		bool writeFileTAB( string sTabFilename );

		bool writeFileXML( string sXMLFilename );

	private:

		bool setupTemplateProteinCombined( string sProRataCombineXMLfilename );

		bool compileLocusList();
		
		string sCombineXMLfilename;
		vector< string > vsLocusList;
		vector< string > vsDescriptionList;
		vector< ProteinCombined * > vpProteinCombined;
		vector< ProteomeInfo * > vpProteomeInfo;

		string sNormalizationPrintOut;

		// the template ProteinCombined to be copied for vpProteinCombined
		ProteinCombined templateProteinCombined;

		// retrieve text from a XML node
		string getValue( TiXmlDocument &txdDoc, const vector<string> &vsTagList );

};
#endif //PROJECTINFO_H
