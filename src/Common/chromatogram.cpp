#include "chromatogram.h"

MS2Scoring::MS2Scoring()
{
	fScore = 0;
	iMSMSscan = 0;
	fRetentionTime = -999;
	sIDfilename = "";
}

Protein::Protein()
{
	sLocus = "";
	sDescription = "";
}

Identification::Identification()
{
	iFirstMS2 = 0;
	iLastMS2 = 0;
	sSequence = "";
	iChargeState = 0;
	
}

Chromatogram::Chromatogram()
{
	iIdentifier = 0;
	sMSfilename = "";
	iScanCount = 0;
	bValidity = true;
	myID.iFirstMS2 = 0;
	myID.iLastMS2 = 0;
	myID.sSequence = "";
	myID.iChargeState = 0;
}

Chromatogram::~Chromatogram()
{

}

int Chromatogram::getIdentifier()
{
	return iIdentifier;
}

const Identification & Chromatogram::getID()
{
	return myID;
}

bool Chromatogram::getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription )
{
	for( int i = 0; i < myID.vProtein.size(); ++i )
	{
		vsLocus.push_back( myID.vProtein[i].sLocus );
		vsDescription.push_back( myID.vProtein[i].sDescription );
	}
	
	return true;

}

float Chromatogram::getMaximumScore()
{
	float fMaximumScore = 0;
	for( int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		if( myID.vMS2Scoring[i].fScore > fMaximumScore )
			fMaximumScore = myID.vMS2Scoring[i].fScore;
	}
	return fMaximumScore;

}

void Chromatogram::getMS2Time( vector< float > & vfMS2Time )
{	
	vfMS2Time.clear();
	for( int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		vfMS2Time.push_back( myID.vMS2Scoring[i].fRetentionTime );
	}
}

vector< unsigned long int > Chromatogram::getMS2ScanNumber()
{	
	vector< unsigned long int > viMS2ScanNumber;
	for( int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		viMS2ScanNumber.push_back( myID.vMS2Scoring[i].iMSMSscan );
	}
	return viMS2ScanNumber;
}
string Chromatogram::getSequence()
{
	return myID.sSequence;
}

int Chromatogram::getChargeState()
{
	return myID.iChargeState;
}

int Chromatogram::getScanCount()
{
	return iScanCount;
}

vector<unsigned long int> Chromatogram::getScanVector()
{
	return viScan;
}

const vector<float> & Chromatogram::getTimeVector()
{
	return vfTime;
}

bool Chromatogram::getIntensityVector( string sName, vector<double> & vdIntensityOutput )
{
	for( int i = 0; i < vSIC.size(); ++i )
	{
		if( vSIC[i].sName == sName )
		{
			vdIntensityOutput = vSIC[i].vdIntensity;
			return true;
		}
	}

	cout << "ERROR: cannot find the chromatogram for the isotopologue " << sName << endl;

	return false;
}

bool Chromatogram::getMZwindows( string sName,  vector< float > & vfLowerMZ, vector< float > & vfUpperMZ  )
{
	vfLowerMZ.clear();
	vfUpperMZ.clear();
	for( int i = 0; i < vSIC.size(); ++i )
	{
		if( vSIC[i].sName == sName )
		{
			vfLowerMZ = vSIC[i].mzWindows.vfLowerMZ;
			vfUpperMZ = vSIC[i].mzWindows.vfUpperMZ;
			return true;
		}
	}

	cout << "ERROR: cannot find the chromatogram for the isotopologue " << sName << endl;

	return false;
}

vector<string> Chromatogram::getAllSICname()
{
	vector< string > vsName;
	for( int i = 0; i < vSIC.size(); ++i )
	{
		vsName.push_back( vSIC[i].sName );
	}

	return vsName;

}

vector<string> Chromatogram::getAllIDfilename()
{
	vector< string > vsName;
	for( int i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		vsName.push_back( myID.vMS2Scoring[i].sIDfilename );
	}

	return vsName;
}

const vector<SIC> & Chromatogram::getAllSIC()
{
	return vSIC;
}

bool Chromatogram::isValid()
{
	return bValidity;
}


unsigned long int Chromatogram::getFullScan4Time( float fTime )
{
	int iTimeSize = vfTime.size();
	for( int i = 0; i < iTimeSize; ++i )
	{
		if( vfTime[i] > fTime )
		{
			return viScan[i];
		}
	}
	return viScan.back();
}

float Chromatogram::getFullScanTime4Time( float fTime )
{
	int iTimeSize = vfTime.size();
	for( int i = 0; i < iTimeSize; ++i )
	{
		if( vfTime[i] > fTime )
		{
			return vfTime[i];
		}
	}
	return vfTime.back();
}

string Chromatogram::getMSfilename()
{
	return sMSfilename;
}

void Chromatogram::setID( const Identification & idInput )
{ 
	myID = idInput; 
}

void Chromatogram::setIdentifier( int iIdentifierInput )
{ 
	iIdentifier = iIdentifierInput; 
}

void Chromatogram::setMSfilename( string sMSfilenameInput )
{ 
	sMSfilename = sMSfilenameInput; 
}

void Chromatogram::setScan( const vector< unsigned long int > & viScanInput )
{ 
	viScan = viScanInput;
	iScanCount = viScan.size();
}

void Chromatogram::setTime( const vector< float > & vfTimeInput )
{ 
	vfTime = vfTimeInput;
}

void Chromatogram::setvSIC( const vector< SIC > & vSICInput )
{ 
	vSIC = vSICInput; 
}

void Chromatogram::setValidity( bool bValidityInput )
{
	bValidity = bValidityInput;
}

bool Chromatogram::writeToXicFile(  FILE * pFile )
{
	ostringstream ossStream;
	int i;
	TiXmlElement * pElementChro = new TiXmlElement( "CHROMATOGRAM" );
	pElementChro->SetAttribute( "MSfile", sMSfilename.c_str() ); 
	pElementChro->SetAttribute( "identifier", iIdentifier ); 
	
	TiXmlElement * pElementIdentification = new TiXmlElement( "IDENTIFICATION" );

	// add SEQUENCE
	TiXmlElement * pElementSequence = new TiXmlElement( "SEQUENCE" );
	TiXmlText * pTextSequence = new TiXmlText( myID.sSequence.c_str() );
	pElementSequence->LinkEndChild( pTextSequence );
	pElementIdentification->LinkEndChild( pElementSequence );
       	
	// add CHARGE_STATE
	TiXmlElement * pElementChargeState = new TiXmlElement( "CHARGE_STATE" );
	ossStream.str("");
	ossStream <<  myID.iChargeState; 
	TiXmlText * pTextChargeState = new TiXmlText( ossStream.str().c_str() );
	pElementChargeState->LinkEndChild( pTextChargeState );
	pElementIdentification->LinkEndChild( pElementChargeState );

	// add MS2SCORING
	for( i = 0; i < myID.vMS2Scoring.size(); ++i )
	{
		TiXmlElement * pElementMS2scoring = new TiXmlElement( "MS2SCORING" );
		
		TiXmlElement * pElementIDFilename = new TiXmlElement( "ID_FILENAME" );
		TiXmlText * pTextIDFilename = new TiXmlText( myID.vMS2Scoring[i].sIDfilename.c_str() );
		pElementIDFilename->LinkEndChild( pTextIDFilename );
		pElementMS2scoring->LinkEndChild( pElementIDFilename );

		TiXmlElement * pElementScanNumber = new TiXmlElement( "SCAN_NUMBER" );
		ossStream.str("");
		ossStream << myID.vMS2Scoring[i].iMSMSscan;
		TiXmlText * pTextScanNumber = new TiXmlText( ossStream.str().c_str() );
		pElementScanNumber->LinkEndChild( pTextScanNumber );
		pElementMS2scoring->LinkEndChild( pElementScanNumber );
		
		TiXmlElement * pElementRT = new TiXmlElement( "RETENTION_TIME" );
		ossStream.str("");
		ossStream << myID.vMS2Scoring[i].fRetentionTime;
		TiXmlText * pTextRT = new TiXmlText( ossStream.str().c_str() );
		pElementRT->LinkEndChild( pTextRT );
		pElementMS2scoring->LinkEndChild( pElementRT );

		
		TiXmlElement * pElementScore = new TiXmlElement( "SCORE" );
		ossStream.str("");
		ossStream << myID.vMS2Scoring[i].fScore;
		TiXmlText * pTextScore = new TiXmlText( ossStream.str().c_str() );
		pElementScore->LinkEndChild( pTextScore );
		pElementMS2scoring->LinkEndChild( pElementScore );
		
		pElementIdentification->LinkEndChild( pElementMS2scoring );

	}
	// add PROTEIN
	for( i = 0; i < myID.vProtein.size(); ++i )
	{
		TiXmlElement * pElementProtein = new TiXmlElement( "PROTEIN" );
		
		TiXmlElement * pElementLocus = new TiXmlElement( "LOCUS" );
		TiXmlText * pTextLocus = new TiXmlText( myID.vProtein[i].sLocus.c_str() );
		pElementLocus->LinkEndChild( pTextLocus );
		pElementProtein->LinkEndChild( pElementLocus );

		TiXmlElement * pElementDescription = new TiXmlElement( "DESCRIPTION" );
		TiXmlText * pTextDescription = new TiXmlText( myID.vProtein[i].sDescription.c_str() );
		pElementDescription->LinkEndChild( pTextDescription );
		pElementProtein->LinkEndChild( pElementDescription );		

		pElementIdentification->LinkEndChild( pElementProtein );
	}
	// add IDENTIFICATION	
	pElementChro->LinkEndChild( pElementIdentification );

	// add SCAN
	TiXmlElement * pElementScan = new TiXmlElement( "SCAN" );
	pElementScan->SetAttribute( "count", viScan.size() ); 
	ossStream.str("");
	if( viScan.size() > 0 )
	{
		for( i = 0; i < viScan.size()-1  ; ++i)
			ossStream << viScan[i] << ",";
		ossStream << viScan[i];
	}
	else
		ossStream << " ";
	TiXmlText * pTextScan = new TiXmlText( ossStream.str().c_str() );  
	pElementScan->LinkEndChild( pTextScan );
	pElementChro->LinkEndChild( pElementScan );	

	// add TIME
	ossStream.str("");
	if( vfTime.size() > 0 )
	{
		for( i = 0; i < vfTime.size()-1  ; ++i)
			ossStream << vfTime[i] << ",";
		ossStream << vfTime[i];
	}
	else
		ossStream << " ";
	TiXmlElement * pElementTime = new TiXmlElement( "TIME" );
	TiXmlText * pTextTime = new TiXmlText( ossStream.str().c_str() );  
	pElementTime->LinkEndChild( pTextTime );
	pElementChro->LinkEndChild( pElementTime );

	ossStream << setprecision(10);
	// add SIC
	int j;
	if( vSIC.size() < 1 )
	{
		// if there is no SIC defined, add an empty one
		TiXmlElement * pElementSIC = new TiXmlElement( "SIC" );
		pElementSIC->SetAttribute( "name", "unknown" ); 
		TiXmlElement * pElementMZwin = new TiXmlElement( "MZ_WINDOW" );
		TiXmlText * pTextMZwin = new TiXmlText( " " );  
		pElementMZwin->LinkEndChild( pTextMZwin );
		pElementSIC->LinkEndChild( pElementMZwin );
		TiXmlElement * pElementIntensity = new TiXmlElement( "INTENSITY" );
		TiXmlText * pTextIntensity = new TiXmlText( " " ); 
		pElementIntensity->LinkEndChild( pTextIntensity );
		pElementSIC->LinkEndChild( pElementIntensity );
		pElementChro->LinkEndChild( pElementSIC );
	}
	else
	{
		for( i = 0; i < vSIC.size(); ++i )
		{
			TiXmlElement * pElementSIC = new TiXmlElement( "SIC" );
			pElementSIC->SetAttribute( "name", vSIC[i].sName.c_str() ); 

			// add MZ_WINDOW

			for( j = 0; j < vSIC[i].mzWindows.vfUpperMZ.size(); ++j )
			{
				TiXmlElement * pElementMZwin = new TiXmlElement( "MZ_WINDOW" );
				ossStream.str("");
				ossStream  << vSIC[i].mzWindows.vfLowerMZ[j] << ", " << vSIC[i].mzWindows.vfUpperMZ[j];
				TiXmlText * pTextMZwin = new TiXmlText( ossStream.str().c_str() );
				pElementMZwin->LinkEndChild( pTextMZwin );
				pElementSIC->LinkEndChild( pElementMZwin );
			}

			// add INTENSITY
			TiXmlElement * pElementIntensity = new TiXmlElement( "INTENSITY" );
			ossStream.str("");
			if( vSIC[i].vdIntensity.size() > 1 )
			{
				for( j = 0; j < vSIC[i].vdIntensity.size()-1  ; ++j )
					ossStream << vSIC[i].vdIntensity[j] << ",";
				ossStream << vSIC[i].vdIntensity[j];
			}
			else
				ossStream << " ";
			TiXmlText * pTextIntensity = new TiXmlText( ossStream.str().c_str() );
			pElementIntensity->LinkEndChild( pTextIntensity );
			pElementSIC->LinkEndChild( pElementIntensity );

			// add SIC
			pElementChro->LinkEndChild( pElementSIC );
		}	
	}
	

	// the first parameter is a File pointer
	// the second parameter is the number of tab indentation
	pElementChro->Print( pFile, 1 );
	fputs( "\n", pFile );

	delete pElementChro;

	return true;
}
	
bool Chromatogram::readXicElement( TiXmlElement * pElementChro )
{
	string sValue;
	istringstream issStream;
	vector< TiXmlElement * > vpElementVector;
	TiXmlElement * pElementTemp;
	int i;
	
	if( !pElementChro )
	{
		cout << "ERROR: unable to read a chromatogram element. " << endl; 
		return false;
	}

	// read the chromatogram identifier
	sValue = pElementChro->Attribute( "identifier" );
	issStream.clear();
	issStream.str( sValue );
	issStream >> iIdentifier;

	// read the MS filename
	sMSfilename = pElementChro->Attribute( "MSfile" );
	
	// push back the element name in the hierarchical order
	// the top level goes first and the leaf node goes las
	vector<string> vsTagList;

	/*
	 * read the identification
	 */
	
	// read the SEQUENCE;
	vsTagList.clear();
	vsTagList.push_back( "IDENTIFICATION" );
	vsTagList.push_back( "SEQUENCE" );
	myID.sSequence = getValue( pElementChro, vsTagList );

	// read the charge state
	vsTagList.pop_back();
	vsTagList.push_back( "CHARGE_STATE" );
	sValue = getValue( pElementChro, vsTagList );
	issStream.clear();
	issStream.str( sValue );
	issStream >> myID.iChargeState;


	// read the MS2SCORING
	vsTagList.pop_back();
	vsTagList.push_back( "MS2SCORING" );

	vpElementVector.clear();
	vpElementVector = getElement( pElementChro, vsTagList );
	for( i = 0; i < vpElementVector.size(); ++i )
	{
		MS2Scoring currentMS2Scoring;
		vsTagList.clear();
		vsTagList.push_back( "ID_FILENAME" );
		currentMS2Scoring.sIDfilename = getValue( vpElementVector[i], vsTagList );

		vsTagList.pop_back();
		vsTagList.push_back( "SCAN_NUMBER" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		issStream >> currentMS2Scoring.iMSMSscan;

		vsTagList.pop_back();
		vsTagList.push_back( "RETENTION_TIME" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		issStream >> currentMS2Scoring.fRetentionTime;

		vsTagList.pop_back();
		vsTagList.push_back( "SCORE" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		issStream >> currentMS2Scoring.fScore;

		myID.vMS2Scoring.push_back( currentMS2Scoring );
		
	}

	// read PROTEIN
	vsTagList.clear();
	vsTagList.push_back( "IDENTIFICATION" );
	vsTagList.push_back( "PROTEIN" );


	vpElementVector.clear();
	vpElementVector = getElement( pElementChro, vsTagList );
	for( i = 0; i < vpElementVector.size(); ++i )
	{
		Protein currentProtein;
		vsTagList.clear();
		vsTagList.push_back( "LOCUS" );
		currentProtein.sLocus = getValue( vpElementVector[i], vsTagList );

		vsTagList.pop_back();
		vsTagList.push_back( "DESCRIPTION" );
		currentProtein.sDescription = getValue( vpElementVector[i], vsTagList );

		myID.vProtein.push_back( currentProtein );
	}	


	// read SCAN
	vsTagList.clear();
	vsTagList.push_back( "SCAN" );
	vpElementVector = getElement( pElementChro, vsTagList );
	sValue = vpElementVector[0]->Attribute( "count" );
	issStream.clear();
	issStream.str( sValue );
	issStream >> iScanCount;

	sValue = getValue( pElementChro, vsTagList );
	for( i = 0; i < sValue.size(); ++i )
	{
		if( sValue[i] == ',' )
			sValue[i] = '\t';
	}
	issStream.clear();
	issStream.str( sValue );
	unsigned long int iScanTemp;
	while( !( issStream.eof() ) )
	{
		issStream >> iScanTemp;
		viScan.push_back( iScanTemp );
	}	

	// read TIME
	vsTagList.clear();
	vsTagList.push_back( "TIME" );
	sValue = getValue( pElementChro, vsTagList );
	for( i = 0; i < sValue.size(); ++i )
	{
		if( sValue[i] == ',' )
			sValue[i] = '\t';
	}
	issStream.clear();
	issStream.str( sValue );
	float fTimeTemp;
	while( !( issStream.eof() ) )
	{
		issStream >> fTimeTemp;
		vfTime.push_back( fTimeTemp );
	}

	// read SIC
	vsTagList.clear();
	vsTagList.push_back( "SIC" );
	vpElementVector.clear();
	vpElementVector = getElement( pElementChro, vsTagList );
	vector< TiXmlElement * > vpElementMZwindow;
	int j;
	int k;
	float fLowerMZtemp;
	float fUpperMZtemp;
	double dIntensityTemp;
	for( i = 0; i < vpElementVector.size(); ++i )
	{
		SIC currentSIC;
		currentSIC.sName = vpElementVector[i]->Attribute( "name" );
		
		vsTagList.clear();
		vsTagList.push_back( "MZ_WINDOW" );
		vpElementMZwindow = getElement( vpElementVector[i], vsTagList );

		vsTagList.clear();
		for( j = 0; j < vpElementMZwindow.size(); ++j )
		{
			sValue = getValue( vpElementMZwindow[j], vsTagList );
			for( k = 0; k < sValue.size(); ++k )
			{
				if( sValue[k] == ',' )
					sValue[k] = '\t';
			}
			issStream.clear();
			issStream.str( sValue );
			issStream >> fLowerMZtemp >> fUpperMZtemp;
			currentSIC.mzWindows.vfLowerMZ.push_back( fLowerMZtemp );
			currentSIC.mzWindows.vfUpperMZ.push_back( fUpperMZtemp );
		}

		vsTagList.clear();
		vsTagList.push_back( "INTENSITY" );
		sValue = getValue( vpElementVector[i], vsTagList );
		for( k = 0; k < sValue.size(); ++k )
		{
			if( sValue[k] == ',' )
				sValue[k] = '\t';
		}

		issStream.clear();
		issStream.str( sValue );
		while( !( issStream.eof() ) )
		{
			issStream >> dIntensityTemp;
			currentSIC.vdIntensity.push_back( dIntensityTemp );
		}
	
		vSIC.push_back( currentSIC );
	}
	
	
	return true;

}

/*
 * this function is copy-and-pasted from proRataConfig.cpp's getValue function
 * the only difference is a TiXmlElement pointer rather than a TiXmlDocument is
 * given as a parameter and the only change is line: txnTemp = pElement->FirstChild( (*itrTagListItr ).c_str() );
 */

/*
 * An utility method to safely extract information from an XML tag.
 * the text inside an XML element is extracted and pasted together if separated
 * by comments or elements.
 * the element is reached by giving a vector of the element name in the order
 * of their hierarchy. Arbitrary number of level can be reached  
 */


string Chromatogram::getValue( TiXmlElement * pElement, const vector<string> &vsTagList )
{

	string sTemp = "";
	if( !pElement )
		return sTemp;
	
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string.

	// creat a pointer to a node
	TiXmlNode * txnTemp = NULL;

	TiXmlText *txs;


	// check if the tree path is not empty
	if ( vsTagList.size() == 0 )
	{
		for( txnTemp = pElement->FirstChild(); txnTemp; txnTemp = txnTemp->NextSibling() )
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
	
	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	/*
	 * the lins is changed from this function's implementation for ProRataConfig
	 */
	// move the pointer to txddoc's first child with the specified tag name
	txnTemp = pElement->FirstChild( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! txnTemp )
	{
		cout << "ERROR: TAG\"" << (*itrTagListItr) << 
			"\" not found in the xic file." << endl;
		return string("");
	}
	itrTagListItr++;
	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		txnTemp = txnTemp->FirstChild( (*itrTagListItr ).c_str() );

		if ( ! txnTemp )
		{
			cout << "ERROR: TAG\"" << (*itrTagListItr) << 
				"\" not found in the xic file." << endl;
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



vector< TiXmlElement * > Chromatogram::getElement(  TiXmlElement * pElement, const vector<string> &vsTagList )
{
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string

	vector< TiXmlElement * > vpElementVector;
	
	if( !pElement )
		return vpElementVector;
	
	// check if the tree path is not empty
	if ( vsTagList.size() < 1 )
		return vpElementVector;
	
	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	/*
	 * the lins is changed from this function's implementation for ProRataConfig
	 */
	// move the pointer to txddoc's first child with the specified tag name
	pElement = pElement->FirstChildElement( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! pElement )
	{
		cout << "ERROR: TAG\"" << (*itrTagListItr) << 
			"\" not found in the xic file." << endl;
		return vpElementVector;
	}

	itrTagListItr++;

	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		pElement = pElement->FirstChildElement( (*itrTagListItr).c_str() );

		if ( ! pElement )
		{
			cout << "ERROR: TAG\"" << (*itrTagListItr) << 
				"\" not found in the xic file." << endl;
			return vpElementVector;
		}

	}

	itrTagListItr--;

	for( pElement = pElement; pElement; pElement = pElement->NextSiblingElement( (*itrTagListItr).c_str() ) )
	{
		vpElementVector.push_back( pElement );
	}

	return vpElementVector;
	


}


bool Chromatogram::readChroFile( const string & sFilename )
{

	// parse the filename to get the MS filename
	sMSfilename = sFilename;
	
	string sLine = "";

	int iScan;
	float fTime;
	double dSample;
	double dReference;
	
	// open the chro file
	ifstream ifsInputFile( sFilename.c_str() );
	if( !ifsInputFile )
	{
		cout << "ERROR: Cannot open the chro file: " << sFilename << endl;
		return false;
	}

	string::size_type startPosition;
	string::size_type endPosition;
	int iLength;
	int i;

	// one chro file contains only one MS2Scoring
	MS2Scoring tempMS2Scoring;
	myID.vMS2Scoring.push_back( tempMS2Scoring );
	myID.vMS2Scoring[0].sIDfilename = sFilename;
	myID.iFirstMS2 = 0;
	myID.iLastMS2 = 0;

	// one chro file have 2 SIC, each of which has one m/z window
	SIC sicSample;
	sicSample.sName = "sample";
	SIC sicRef;
	sicRef.sName = "reference";
	while( !ifsInputFile.eof() )
	{
		getline(ifsInputFile, sLine);

		if ( sLine.find( "Locus:", 0 ) == 0 )
		{
			
			startPosition = 6;
			endPosition = sLine.find( ",", startPosition ); 
			while( endPosition != string::npos )
			{
				Protein currentProtein;
				// one character shorter to remove ","
				iLength = endPosition - startPosition;
				istringstream issContent( sLine.substr( startPosition, iLength ) );
				issContent >> currentProtein.sLocus;
				myID.vProtein.push_back( currentProtein );
				endPosition++;
				startPosition = ( endPosition + 1 );
				endPosition = sLine.find( ",", startPosition );
			}
			Protein currentProtein;
			iLength = sLine.length() - startPosition;
			istringstream issContent( sLine.substr( startPosition, iLength ) );
			issContent >> currentProtein.sLocus;
			myID.vProtein.push_back( currentProtein );
			continue;
		}
		
		if ( sLine.find( "Description:", 0 ) == 0 )
		{
			// if there is more than one locus, there is no way to parse out individual locus'
			// description with comma from this line. Just leave the description for those locuses
			// empty
			if( ! ( myID.vProtein.size() == 1 ) )
				continue;
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			issContent >> myID.vProtein[0].sDescription;
			continue;
		}

		

		if ( sLine.find( "MS/MS Scan Number:", 0 ) == 0 )
		{
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			issContent >> myID.vMS2Scoring[0].iMSMSscan;
			continue;
		}
		
		if ( sLine.find( "XCorr:", 0 ) == 0 )
		{
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			issContent >> myID.vMS2Scoring[0].fScore;
			continue;
		}
		
		if ( sLine.find( "Sequence:", 0 ) == 0 )
		{
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			issContent >> myID.sSequence;
			continue;	
		}

		if ( sLine.find( "Charge:", 0 ) == 0 )
		{
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			issContent >> myID.iChargeState;
			continue;	
		}

		if ( sLine.find( "M/Z Range Sample:", 0 ) == 0 )
		{
			for( i = 0; i < sLine.length(); ++i )
			{
				if( sLine[i] == '-' )
					sLine[i] = '\t';
			}
					
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			float fLowerMZtemp;
			float fUpperMZtemp;
			issContent >> fLowerMZtemp >> fUpperMZtemp;
			sicSample.mzWindows.vfUpperMZ.push_back( fUpperMZtemp );
			sicSample.mzWindows.vfLowerMZ.push_back( fLowerMZtemp );
			continue;	
		}
		
		if ( sLine.find(  "M/Z Range Reference:", 0 ) == 0 )
		{
			for( i = 0; i < sLine.length(); ++i )
			{
				if( sLine[i] == '-' )
					sLine[i] = '\t';
			}
					
			istringstream issContent(sLine);
			issContent.seekg( sLine.find( ":" ) + 1 );
			float fLowerMZtemp;
			float fUpperMZtemp;
			issContent >> fLowerMZtemp >> fUpperMZtemp;
			sicRef.mzWindows.vfUpperMZ.push_back( fUpperMZtemp );
			sicRef.mzWindows.vfLowerMZ.push_back( fLowerMZtemp );
			continue;	
		}
		
		// Chromatogram data
		if ( sLine.compare( 0, 15, "[CHROMATOGRAMS]" ) == 0 )
		{
			getline(ifsInputFile, sLine);
			while( !  ifsInputFile.eof()  )
			{
				ifsInputFile >> iScan;

				viScan.push_back( iScan );

				ifsInputFile >> fTime;
				vfTime.push_back( fTime );

				ifsInputFile >> dSample;
				sicSample.vdIntensity.push_back( dSample );

				ifsInputFile >> dReference;
				sicRef.vdIntensity.push_back( dReference );

			}
		} // end Chromatogram data
	}
	
	vSIC.push_back( sicRef );
	vSIC.push_back( sicSample );
	iScanCount = viScan.size();
	setValidity( true );

	unsigned long int iScanDifference = abs( (long)(viScan[0] - myID.vMS2Scoring[0].iMSMSscan) );
	unsigned long int iMinScanDifference = iScanDifference;
	int index = 0;
	for( i = 1; i < viScan.size(); ++i )
	{
		iScanDifference = abs( (long)( viScan[i] - myID.vMS2Scoring[0].iMSMSscan ) );
		if( iScanDifference < iMinScanDifference )
		{
			index = i;
			iMinScanDifference = iScanDifference;
		}
	}

	myID.vMS2Scoring[0].fRetentionTime = vfTime[index];
	
	return true;

}


