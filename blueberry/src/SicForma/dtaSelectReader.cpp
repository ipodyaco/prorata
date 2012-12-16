#include "dtaSelectReader.h"

DTASelectReader::DTASelectReader()
{
	// constructor
}

DTASelectReader::~DTASelectReader()
{
	// destructor
}

list< Identification * > DTASelectReader::getIDlist( string sFilename )
{
	
	list< Identification * > lpIDlist;
	
	getIDlist( sFilename, lpIDlist );

	return lpIDlist;

}

bool DTASelectReader::getIDlist( string sFilename, list< Identification * > & lpIDlist )
{
	string sCurrentLine;

	map<string, Identification *> mapScan2ID;
	map<string, Identification *>::iterator mapIterScan2ID;
	
	// intermediate preceding protein line
	Protein precedingProtein;
	
	// open the input DTASelect-filter file
	ifstream fileDTASelect( sFilename.c_str() );
	if( !fileDTASelect )
	{
		cout << "ERROR: Cannot open the .pro2psm.txt file: " << sFilename << endl;
		return false;
	}

	unsigned int iFieldProteinID = 0;
	unsigned int iFieldProteinDescription = 0;
	// move to the protein header lines, which is the first + line
	while( getline(fileDTASelect, sCurrentLine) )
	{
		if( sCurrentLine[0] == '+' )
		{
			TokenVector words(sCurrentLine, "\t");
			if(!findField(words, "ProteinID", iFieldProteinID))
				return false;
			if(!findField(words, "ProteinDescription", iFieldProteinDescription))
				return false;
			break;
		}
	}

	unsigned int iFieldFilename = 0;
	unsigned int iFieldScanNumber = 0;
	unsigned int iFieldParentCharge = 0;
	unsigned int iFieldSearchName = 0;
	unsigned int iFieldScore = 0;
	unsigned int iFieldIdentifiedPeptide = 0;
	// move to the peptide header lines, which is the first * line
	while( getline(fileDTASelect, sCurrentLine) )
	{
		if( sCurrentLine[0] == '*' )
		{
			TokenVector words(sCurrentLine, "\t");
			if(!findField(words, "Filename", iFieldFilename))	return false;
			if(!findField(words, "ScanNumber", iFieldScanNumber))	return false;
			if(!findField(words, "ParentCharge", iFieldParentCharge))	return false;
			if(!findField(words, "SearchName", iFieldSearchName))	return false;
			if(!findField(words, "Score", iFieldScore))		return false;
			if(!findField(words, "IdentifiedPeptide", iFieldIdentifiedPeptide))	return false;
			break;
		}
	}
	


	// move down the file line by line
	while( getline(fileDTASelect, sCurrentLine) )
	{
		if( sCurrentLine[0] == '+' )
		{
			TokenVector words(sCurrentLine, "\t");			
			precedingProtein.sLocus = words[iFieldProteinID];
			precedingProtein.sDescription = words[iFieldProteinDescription];
		}
		else if( sCurrentLine[0] == '*' )
		{
			TokenVector words(sCurrentLine, "\t");
			string sFilename = words[iFieldFilename];
			string sScanNumber = words[iFieldScanNumber];
			string sParentCharge = words[iFieldParentCharge];
			string sSearchName = words[iFieldSearchName];
			string sScore = words[iFieldScore];
			string sIdentifiedPeptide = words[iFieldIdentifiedPeptide];

			string sCurrentScan = sFilename + "." + sScanNumber + "." + 
				sParentCharge + "." + sSearchName + "." + sIdentifiedPeptide;

			mapIterScan2ID = mapScan2ID.find(sCurrentScan);
			if(mapIterScan2ID != mapScan2ID.end())
			{
				// this scan is non-unique and has been added
				// just add this protein to this scan too
				mapIterScan2ID->second->vProtein.push_back(precedingProtein);

			}
			else
			{
				// creat a Identification instance for this peptide line
				Identification* pID = new Identification;
				pID->sSequence = sIdentifiedPeptide;
				pID->iChargeState = atoi( sParentCharge.c_str() );
				MS2Scoring currentMS2Scoring;
				currentMS2Scoring.sIDfilename = sSearchName + "." + sFilename + "." + sScanNumber + "." + sParentCharge ;
				currentMS2Scoring.fScore = ( float )atof( sScore.c_str() );
				currentMS2Scoring.iMSMSscan = atoi( sScanNumber.c_str() );
				pID->vMS2Scoring.push_back( currentMS2Scoring );
				pID->vProtein.push_back(precedingProtein);
				mapScan2ID[sCurrentScan] = pID;
			}

		}
		else
		{
			// ignore unknown line
		}

	}


	for (mapIterScan2ID = mapScan2ID.begin(); mapIterScan2ID != mapScan2ID.end(); ++mapIterScan2ID)
	{
		lpIDlist.push_back(mapIterScan2ID->second);
	}

	// if nothing is saved to lpIDlist, return false
	if( lpIDlist.size() > 0 )
		return true;
	else
		return false;

}

bool DTASelectReader::findField(const TokenVector & words, string targetWord, unsigned int & iFieldNumber)
{
	iFieldNumber = 0;
	unsigned int i = 0;
	for( i = 0; i < words.size(); i++ )
	{
		if(words[i] == targetWord)
		{
			iFieldNumber = i;
			return true;
		}
	}
	cout << "Error: Cannot find the field: " << targetWord << endl;
	return false;
}

