
#include "proRataConfig.h"

string ProRataConfig::sFilename = "ProRataConfig.xml";

bool ProRataConfig::bIfWriteChro = false;
bool ProRataConfig::bIsLabelFree = false;

#if _WIN32
	string ProRataConfig::sWorkingDirectory = ".\\";
#else
	string ProRataConfig::sWorkingDirectory = ".";
#endif

// variables from the SIC_EXTRACTION element
float ProRataConfig::fMinutesBeforeMS2 = 2;
float ProRataConfig::fMinutesAfterMS2 = 2;
float ProRataConfig::fMinutesBetweenMS2 = 0.5;
double ProRataConfig::dPlusMZerror = 0.5;
double ProRataConfig::dMinusMZerror = 0.5;
double ProRataConfig::dIsotopicEnvelopCutoff = 0.1;

// variables from the PEPTIDE_QUANTIFICATION element
int ProRataConfig::iPeakDetectionFOrder = 2;
int ProRataConfig::iPeakDetectionWindowSize = 7;

float ProRataConfig::fLeftPeakShift = 0;
float ProRataConfig::fRightPeakShift = 0;

string ProRataConfig::sNumeratorIsotopologue = "treatment";
string ProRataConfig::sDenominatorIsotopologue = "reference";

double ProRataConfig::dPCAMinLog2Ratio = -7;
double ProRataConfig::dPCAMaxLog2Ratio = 7;

// variables from the PROTEIN_QUANTIFICATION element

bool ProRataConfig::bRemoveAmbiguousPeptide = true;

int ProRataConfig::iMinPeptideNumber = 2;
double ProRataConfig::dMaxCIwidth = 7;

double ProRataConfig::dMinLog2SNR = 1;
double ProRataConfig::dMaxLog2SNR = 4;

double ProRataConfig::dMLEMinLog2Ratio = -7;
double ProRataConfig::dMLEMaxLog2Ratio = 7;

double ProRataConfig::dLog2RatioDiscretization = 0.1;

double ProRataConfig::dSDSlope = -0.26;
double ProRataConfig::dSDIntercept = 1.3;

double ProRataConfig::dMeanSlope = 1.2;
double ProRataConfig::dMeanIntercept = 0;

double ProRataConfig::dSmoothingProbSpace = 0.15;
double ProRataConfig::dLnLikelihoodCutoffOffset = 1.96;

ProRataConfig* ProRataConfig::proRataConfigSingleton = 0;

ProRataConfig::ProRataConfig()
{}

bool ProRataConfig::setFilename( const string & sConfigFileName )
{
	if ( proRataConfigSingleton == 0 )
	{
		proRataConfigSingleton = new ProRataConfig;
	}

	sFilename = sConfigFileName;
	
	// Creat a TinyXML document for ProRataConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	proRataConfigSingleton->getParameters( txdConfigFile );
	proRataConfigSingleton->setIsLabelFree();	
	// If everything goes fine return 0.
	return true;

}

bool ProRataConfig::setWorkingDirectory( const string & sDirectoryName )
{
	if( sDirectoryName[ sDirectoryName.size() - 1 ] == ProRataConfig::getSeparator() )
	{
		sWorkingDirectory = sDirectoryName;
	}
	else
	{
		sWorkingDirectory = sDirectoryName + ProRataConfig::getSeparator();
	}

	return true;
}	

char ProRataConfig::getSeparator()
{
#if _WIN32
		return '\\' ;
#else
		return '/' ;
#endif
}



bool ProRataConfig::getAtomIsotopicComposition( char cAtom, 
		vector<double> & vdAtomicMass,  
		vector<double> & vdNaturalComposition,
		vector<double> & vdEnrichedComposition)
{

	// clear the input vectors
	vdAtomicMass.clear();
	vdNaturalComposition.clear();
	vdEnrichedComposition.clear();

	// Creat a TinyXML document for ProRataConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	string sData;
	istringstream issStream;
	double dValue;
	string sAtom = "X";
	sAtom[0] = cAtom;
	
	// creat the path to the MASS_DA element of the input sAtom
	vector<string> vsTagList;
	vsTagList.push_back( "CONFIG" );
	vsTagList.push_back( "SIC_EXTRACTION" );
	vsTagList.push_back( "ATOM_ISOTOPIC_COMPOSITION" );
	vsTagList.push_back( sAtom );
	vsTagList.push_back( "MASS_DA" );

	// get the text inside and extract the value
	sData = getValue( txdConfigFile, vsTagList );
	replaceDelimitor( sData, ',', '\t' );
	// clear end of file state 
	issStream.clear();
	// re-set the string associated with issStream
	issStream.str( sData );
	while( !( issStream.eof() ) )
	{
		issStream >> dValue;
		vdAtomicMass.push_back( dValue );
	}

	// move to the NATURAL element
	vsTagList.pop_back();
	vsTagList.push_back( "NATURAL" );
	sData = getValue( txdConfigFile, vsTagList );
	replaceDelimitor( sData, ',', '\t' );
	issStream.clear();
	issStream.str( sData );
	while( !( issStream.eof() ) )
	{
		issStream >> dValue;
		vdNaturalComposition.push_back( dValue );
	}
	
	// move to the ENRICHED element
	vsTagList.pop_back();
	vsTagList.push_back( "ENRICHED" );
	sData = getValue( txdConfigFile, vsTagList );
	replaceDelimitor( sData, ',', '\t' );
	issStream.clear();
	issStream.str( sData );
	while( !( issStream.eof() ) )
	{
		issStream >> dValue;
		vdEnrichedComposition.push_back( dValue );
	}

	return true;
}

bool ProRataConfig::getResidueAtomicComposition(residueMap & mIsotopologue)
{
	mIsotopologue.clear();

	// Creat a TinyXML document for ProRataConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}	

	/*
	 * move the node pointer to RESIDUE_ATOMIC_COMPOSITION
	 * if failed, return an empty map
	 */
	TiXmlNode * txnTemp = NULL;

	txnTemp = txdConfigFile.FirstChild( "CONFIG" );
	if ( ! txnTemp )
		return false;
	
	txnTemp = txnTemp->FirstChild( "SIC_EXTRACTION" );
	if ( ! txnTemp )
		return false;
	
	txnTemp = txnTemp->FirstChild( "RESIDUE_ATOMIC_COMPOSITION" );
	if ( ! txnTemp )
		return false;

	// a node pointer to a ISOTOPOLOGUE element
	TiXmlNode * txnTable = NULL;
	
	// a element pointer for retrieving the attribute of the ISOTOPOLOGUE element
	TiXmlElement * txeElement = NULL;

	// a text pointer for retrieving the text inside ISOTOPOLOGUE element
	TiXmlText * txsText = NULL;
		
	/*
	 * the node type can any one of the following:
	 * enum NodeType
	 * {
	 * 	DOCUMENT,
	 * 	ELEMENT,
	 * 	COMMENT,
	 * 	UNKNOWN,
	 * 	TEXT,
	 * 	DECLARATION,
	 * 	TYPECOUNT
	 * };
	 * 
	*/
	
	// loop thru all ISOTOPOLOGUE elements	
	for( txnTemp = txnTemp->FirstChild( "ISOTOPOLOGUE" ); txnTemp; txnTemp = txnTemp->NextSibling( "ISOTOPOLOGUE" )  )
	{
		// the string for holding the text inside
		string sTable = "";

		/*
		 * loop thru all the text nodes inside ISOTOPOLOGUE;
		 * the text nodes can be separated by the comment nodes or other
		 * cast the node to a text node, only if it is of type TEXT, which equals 4
		 * then concatenate all the text
		 */
		for( txnTable = txnTemp->FirstChild("R"); txnTable; txnTable = txnTable->NextSibling("R") )
		{
			txsText =  txnTable->FirstChild()->ToText();
			string sResidueLine = txsText->Value();

			string symbolName;
			string symbol;
			size_t position;

			// if the XML reserved characters, & and <, are used to represent PTMs in the peptide,
			// they can specified as amp and lt in the config file
			// here they will be changed back to symbols

			// replace amp with &
			symbolName = "amp";
			symbol = "&";
			position = sResidueLine.find(symbolName);
			if( position != std::string::npos )
				sResidueLine.replace(position, symbolName.length(), symbol);

			// replace lt with <
			symbolName = "lt";
			symbol = "<";
			position = sResidueLine.find(symbolName);
			if( position != std::string::npos )
				sResidueLine.replace(position, symbolName.length(), symbol);

			sTable.append( sResidueLine );
			sTable.append( "\n" );
		}
		replaceDelimitor( sTable, ',', '\t' );

//		cout << "Config Residue Table" << endl;
//		cout << sTable << endl;

		// points txeElement to the ISOTOPOLOGUE element 
		txeElement = txnTemp->ToElement();
		if( txeElement )
		{
			// get the "name" Attribute of the element
			string sName = txeElement->Attribute( "name" );
			// save the sName-sTable pair to the return value
			mIsotopologue[ sName ] = sTable;
		}
		
	}

	
			
	return true;

}

void ProRataConfig::setIsLabelFree()
{
	residueMap mIsotopologue;
	if(!getResidueAtomicComposition(mIsotopologue)){
		bIsLabelFree = false;
	}

	if(mIsotopologue.size() == 1){
		bIsLabelFree = true;
	}
	else{
		bIsLabelFree = false;
	}
}


void ProRataConfig::getParameters( TiXmlDocument & txdConfigFile )
{

	// strings used to specify the path
	string sMainTag = "CONFIG";
	string sModuleTag;

	// push back the element name in the hierarchical order
	// the top level goes first and the leaf node goes last
	vector<string> vsTagList;	
	
	string sTemp;
	istringstream issStream;

	// Extract the elements inside <SIC_EXTRACTION>
	sModuleTag = "SIC_EXTRACTION";
	
	vsTagList.clear();
	vsTagList.push_back( sMainTag );
	vsTagList.push_back( sModuleTag );
	vsTagList.push_back( "" );
	
	vsTagList[2] = "RETENTION_TIME_INTERVAL";
	vsTagList.push_back( "MINUTES_BEFORE_MS2" );
	sTemp = getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> fMinutesBeforeMS2;
	
	vsTagList[3] = "MINUTES_AFTER_MS2";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> fMinutesAfterMS2;

	vsTagList[3] = "MINUTES_BETWEEN_DUPLICATE_MS2";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> fMinutesBetweenMS2 ;
	
	vsTagList[2] = "MASS_TO_CHARGE_INTERVAL";
	vsTagList[3] = "PLUS_MZ_ERROR";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dPlusMZerror;

	vsTagList[3] = "MINUS_MZ_ERROR";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMinusMZerror;

	vsTagList[3] = "ISOTOPIC_ENVELOP_CUTOFF";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dIsotopicEnvelopCutoff;
	
	// extract elements inside PEPTIDE_QUANTIFICATION
	sModuleTag = "PEPTIDE_QUANTIFICATION";
	vsTagList.clear();
	vsTagList.push_back( sMainTag );
	vsTagList.push_back( sModuleTag );
	vsTagList.push_back( "PEAK_DETECTION" );
	vsTagList.push_back( "CHROMATOGRAM_SMOOTHING" );
	vsTagList.push_back( "WINDOW_SIZE" );
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iPeakDetectionWindowSize;
	
	vsTagList[4] = "ORDER";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iPeakDetectionFOrder;
	
	vsTagList[3] = "PEAK_SHIFT";
	vsTagList[4] = "LEFT";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> fLeftPeakShift;

	vsTagList[4] = "RIGHT";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> fRightPeakShift;
	
	vsTagList.pop_back();
	vsTagList[2] = "LOG2_RATIO";
	vsTagList[3] = "MINIMUM";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dPCAMinLog2Ratio;

	vsTagList[3] = "MAXIMUM";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dPCAMaxLog2Ratio;

	vsTagList[2] = "ABUNDANCE_RATIO";
	vsTagList[3] = "NUMERATOR_ISOTOPOLOGUE";
	sNumeratorIsotopologue = getValue( txdConfigFile, vsTagList ) ;

	vsTagList[3] = "DENOMINATOR_ISOTOPOLOGUE";
	sDenominatorIsotopologue = getValue( txdConfigFile, vsTagList ) ;

	vsTagList.pop_back();
	vsTagList[2] = "LOG2_SNR_CUTOFF";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMinLog2SNR;	
	
	vsTagList[2] = "REMOVE_AMBIGUOUS_PEPTIDES";
	sTemp = getValue( txdConfigFile, vsTagList );
	if( sTemp == "true" || sTemp == "True" || sTemp == "TRUE" || sTemp == "T" )
		bRemoveAmbiguousPeptide = true;
	else
		bRemoveAmbiguousPeptide = false;
	
	// extract elements inside PROTEIN_QUANTIFICATION
	sModuleTag = "PROTEIN_QUANTIFICATION";
	vsTagList.clear();
	vsTagList.push_back( sMainTag );
	vsTagList.push_back( sModuleTag );

	vsTagList.push_back( "MIN_PEPTIDE_NUMBER" );
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iMinPeptideNumber;

	vsTagList[2] = "LOG2_RATIO_DISCRETIZATION";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dLog2RatioDiscretization;

	vsTagList[2] = "SMOOTHING_PROBABILITY_SPACE";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dSmoothingProbSpace;

	vsTagList[2] = "LN_LIKELIHOOD_CUTOFF_OFFSET";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dLnLikelihoodCutoffOffset;
	
	vsTagList[2] = "MAX_CI_WIDTH";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMaxCIwidth;
	
	vsTagList[2] = "MAX_LOG2_SNR";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMaxLog2SNR;	
	
	vsTagList[2] = "LOG2_RATIO";
	vsTagList.push_back( "MINIMUM" );
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMLEMinLog2Ratio;	

	vsTagList[3] = "MAXIMUM";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMLEMaxLog2Ratio;	

	vsTagList[2] = "STANDARD_DEVIATION";
	vsTagList[3] = "SLOPE";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dSDSlope;

	vsTagList[3] = "INTERCEPT";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dSDIntercept;

	vsTagList[2] = "MEAN";
	vsTagList[3] = "SLOPE";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMeanSlope;

	vsTagList[3] = "INTERCEPT";
	sTemp =  getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMeanIntercept;

}


/*
 * An utility method to safely extract information from an XML tag.
 * the text inside an XML element is extracted and pasted together if separated
 * by comments or elements.
 * the element is reached by giving a vector of the element name in the order
 * of their hierarchy. Arbitrary number of level can be reached  
 */

string ProRataConfig::getValue( TiXmlDocument &txdDoc, 
		const vector<string> &vsTagList )
{
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string.

	// creat a pointer to a node
	TiXmlNode * txnTemp = NULL;

	// check if the tree path is not empty
	if ( vsTagList.size() < 1 )
		return string("");
	
	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	// move the pointer to txddoc's first child with the specified tag name
	txnTemp = txdDoc.FirstChild( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! txnTemp )
	{
	//	cout << "ERROR: TAG \"" << (*itrTagListItr) << "\" not found in the configuration file." << endl;
		return string("");
	}

	itrTagListItr++;

	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		txnTemp = txnTemp->FirstChild( (*itrTagListItr ).c_str() );

		if ( ! txnTemp )
		{
		//	cout << "ERROR: TAG \"" << (*itrTagListItr) << "\" not found in the configuration file." << endl;
			return string("");
		}

	}

	// move the iterator back to point it to the last element name
	itrTagListItr--;

	/*
	 * inside the pointed element, there could be a mixture of
	 * text nodes, comment nodes and element nodes
	 * loop thru every nodes and for each text nodes, retrieve their text
	 * concatenate the text together and return them
	 */
	
	TiXmlText *txs;
	string sTemp = "";
	
	// point txnTemp to the child nodes and loop thru every child node
	for( txnTemp = txnTemp->FirstChild(); txnTemp; txnTemp = txnTemp->NextSibling() )
	{
		// if this node is pointing to a node of type TEXT, which equals 4 in enum NodeType
		if( txnTemp->Type() == 4 )
		{
			// cast txnTemp to a text node
			txs = txnTemp->ToText();
			// get txnTemp's value and then append it to sTemp
			if( txs )
				sTemp.append( txs->Value() );
		}
	}

	return sTemp;
}

void ProRataConfig::replaceDelimitor( string & sLine, char cOldDelimitor, char cNewDelimitor )
{
	int iLength = sLine.length();
	for( int i = 0; i < iLength; ++i )
	{
		if( sLine[i] == cOldDelimitor )
			sLine[i] = cNewDelimitor;
	}
	return;
}

bool ProRataConfig::writeConfigXML( string sOutputFile, int iTabDepth, bool bSicExtract, bool bPeptideQuan, bool bProteinQuan )
{
	FILE * pFILE;
	if( ( pFILE = fopen( sOutputFile.c_str(), "a" ) ) == NULL  ) 
	{
		cout << "ERROR: cannot write Config XML to the file: " << sOutputFile << endl;
		return false;
	}

//	TiXmlElement * pElementConfig = new TiXmlElement( "CONFIG" );
//	pElementConfig->SetAttribute( "version", getProRataVersion() );
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR: Loading Configuration file!" << endl;
		return false;
	}

	TiXmlElement * pElementConfig = txdConfigFile.FirstChildElement( "CONFIG" );

	if( !pElementConfig )
	{
		cout <<  "ERROR: cannot find the CONFIG element in the given Configuration file! " << endl;
		return false;
	}
	
	TiXmlNode * pElementSIC = pElementConfig->FirstChild( "SIC_EXTRACTION" );
	TiXmlNode * pElementPeptideQuan = pElementConfig->FirstChild( "PEPTIDE_QUANTIFICATION" );
       	TiXmlNode * pElementProteinQuan = pElementConfig->FirstChild( "PROTEIN_QUANTIFICATION" );

	
	if( (!bSicExtract) && pElementSIC )
	{
		pElementConfig->RemoveChild(pElementSIC);
	}
	
	if( (!bPeptideQuan) && pElementPeptideQuan )
	{
		pElementConfig->RemoveChild(pElementPeptideQuan);
	}

	if( (!bProteinQuan) && pElementProteinQuan )
	{
		pElementConfig->RemoveChild(pElementProteinQuan);
	}
	
	pElementConfig->Print( pFILE, iTabDepth );
	fputs( "\n", pFILE );

	fclose( pFILE );
	return true;
}
