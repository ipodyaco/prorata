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
	
	// a flag for whether the preceding line is a protein line 
	bool bIsLastLineProteinLine = false;

	// a protein vector containing the parsed proteins
	vector< Protein > vProtein;
	
	// open the input DTASelect-filter file
	ifstream fileDTASelect( sFilename.c_str() );
	if( !fileDTASelect )
	{
		cout << "ERROR: Cannot open the DTASelect-filter file: " << sFilename << endl;
		return false;
	}

	// move to the start line
	while( getline(fileDTASelect, sCurrentLine) )
	{
		if( whatIsThisLine( sCurrentLine ) == startLine )
			break;
	}

	// move down the file line by line
	while( getline(fileDTASelect, sCurrentLine) )
	{
		// check if this is the last data line
		if( whatIsThisLine( sCurrentLine ) == endLine )
			break;
		
		// determine what is this line
		if( whatIsThisLine( sCurrentLine ) == proteinLine )
		{
			// if the last line is not a protein line
			// then the protein vector is cleared
			// otherwise the two proteins are grouped into the vector
			if( !bIsLastLineProteinLine )
				vProtein.clear();
			
			// attempt to process the protein line and save to vProtein
			if( processProteinLine( sCurrentLine, vProtein ) )
			{
				// this line is a protein line, which is the last line relative to the next line 
				bIsLastLineProteinLine = true;
			}
			else
			{
				bIsLastLineProteinLine = false;
				cout << "WARNING: the following protein line in DTASelect-filter cannot be parsed and is skipped: " << endl;
				cout << sCurrentLine << endl; 
			}
		}
		else if( whatIsThisLine(sCurrentLine) == peptideLine )
		{
			// creat a Identification instance for this peptide line
			// if this peptide is processed successfully, save it to the vpIDlist
			// otherwise delete it and free the memory
			Identification* pID = new Identification;
			if( processPeptideLine( sCurrentLine, vProtein, pID ) )
			{
				lpIDlist.push_back( pID );
			}
			else
			{
				delete pID;
				cout << "WARNING: the following peptide line in DTASelect-filter cannot be parsed and is skipped: " << endl;
				cout << sCurrentLine << endl;
			}
			// this line is not a protein line, which is the last line relative to the next line 
			bIsLastLineProteinLine = false;
		}
		else
		{
			// cannot determine what is this line
			cout << "WARNING: the following unknown line in DTASelect-filter cannot be parsed and is skipped: " << endl;
			cout << sCurrentLine << endl;
			bIsLastLineProteinLine = false;
		}
		

	}

	// if nothing is saved to lpIDlist, return false
	if( lpIDlist.size() > 0 )
		return true;
	else
		return false;

}

bool DTASelectReader::processProteinLine( string sLine, vector< Protein > & vProtein )
{
	
	Protein currentProtein;

	// move position to the first tab
	string::size_type position = sLine.find( "\t", 0 ) ;
	
	if( position != string::npos )
	{
		// save the characters from the first character to the character proceding the tab
		// to currentProtein.sLocus
		currentProtein.sLocus = sLine.substr( 0, position );
	}
	else
		return false;

	// move position passing seven tabs
	// the protein description is the nineth tab-delimited field in 
	for( int i = 0 ; i < 7; i++ )
	{
		// move position next to the character trailing the found tab
		++position;
		// find the next tab
		position = sLine.find( "\t", position );
		if( position ==  string::npos )
			return false;
	}

	// move position to the first character of the description
	// extract all characters from here to the end
	++position;
	currentProtein.sDescription = sLine.substr( position, ( sLine.length() - position ) );

	// save to vProtein
	vProtein.push_back( currentProtein );
	return true;

}


bool DTASelectReader::processPeptideLine( string sLine, const vector< Protein > & vProtein, Identification* pID )
{
	// member variables for *pID
	MS2Scoring currentMS2Scoring;
	string sSequence;
	int iChargeState;
	
	// the temporary position that is moved step by step	
	string::size_type tempPosition;

	// the start position of a field
	string::size_type startPosition;

	// the end position of a field
	string::size_type endPosition;

	// the length of a field
	int iLength;

	/*
	 *  parse out ID filename
	 *  sIDfilename is the field following the first tab
	 */
	string sIDfilename;

	// find the first tab
	tempPosition = sLine.find( "\t", 0 ) ;
	if( tempPosition  == string::npos )
		return false;

	// the start position of sIDfilename is the one following the first tab
	startPosition = ++tempPosition;

	// move tempPosition to the next tab
	tempPosition = sLine.find( "\t", tempPosition );
	if( tempPosition == string::npos)
		return false;
	
	// calculate the length of sIDfilename
	iLength = tempPosition - startPosition;
	sIDfilename =  sLine.substr( startPosition, iLength );
	currentMS2Scoring.sIDfilename = sIDfilename;

	/*
	 * parse out Xcorr
	 * Xcorr is the next tab-delimited field
	 */
	startPosition = ++tempPosition;
	tempPosition = sLine.find( "\t", tempPosition );
	if( tempPosition == string::npos )
		return false;
	iLength = tempPosition - startPosition;
	string sXcorr = sLine.substr( startPosition, iLength );
	currentMS2Scoring.fScore = ( float )atof( sXcorr.c_str() );

	/*
	 * parse out peptide sequence
	 * the peptide sequence is the last tab-delimited field
	 */
	// find the last tab
	tempPosition = sLine.rfind( "\t", ( sLine.length() - 1 ) );
	startPosition = ++tempPosition;
	// calculate the length
	iLength = sLine.length() - startPosition;
	sSequence = sLine.substr( startPosition, iLength );

	/* parse out charge state from the sIDfilename
	 * the charge state is the last '.'-delimited field
	 */
	string sChargeStateString;
	tempPosition = sIDfilename.rfind( ".", (sLine.length() - 1) );
	if( tempPosition == string::npos )
		return false;
	startPosition = ( tempPosition + 1 );
	iLength = sLine.length() - startPosition;
	sChargeStateString = sIDfilename.substr( startPosition, iLength );
	iChargeState = atoi( sChargeStateString.c_str() );

	/*
	 * parse out scan number from the sIDfilename
	 * the scan number is the next '.'-delimited field preceding the charge state field
	 * this scan number is actually the last scan number, which is
	 * usually the same as the first scan number
	 */
	endPosition = --tempPosition;
	tempPosition = sIDfilename.rfind( ".", tempPosition );
	if( tempPosition == string::npos )
		return false;
	startPosition = ( tempPosition + 1 );
	iLength = endPosition - startPosition + 1;
	string sScanString = sIDfilename.substr( startPosition, iLength );
	currentMS2Scoring.iMSMSscan = atoi( sScanString.c_str() );
	if( currentMS2Scoring.iMSMSscan < 1 )
		return false;

	// save everything to *pID
	// the vProtein is one of the input parameters
	pID->sSequence = sSequence;
	pID->iChargeState = iChargeState;
	pID->vMS2Scoring.push_back( currentMS2Scoring );
	pID->vProtein = vProtein;

	return true;
	
}



LineType DTASelectReader::whatIsThisLine( string sLine )
{
	// if the character is a '*', which mean a unique peptide
	// then this is a peptide line
	if( sLine[0] == '*' )
		return peptideLine;

	// if "Unique" is found as the first word
	// then this a start line
	if( sLine.find( "Unique", 0 ) == 0 )
		return startLine;

	// if "Protein" is found starting from the second position
	// the first position in a ending line is a tab 
	if( sLine.find( "Proteins", 0 ) == 1 )
		return endLine;
			
	/*
	 * a peptide line should has a ID filename, as the field following
	 * to the first tab
	 */
	bool bIsPeptideLine = false;
	// find the first tab
	string::size_type tempPosition;	
	tempPosition = sLine.find( "\t", 0 ) ;
	if( tempPosition  != string::npos )
	{
		// find the second tab
		++tempPosition;
		tempPosition = sLine.find( "\t", tempPosition );
		if( tempPosition != string::npos)
		{
			// a peptide line should has an ID filename, which use '.' to delimit the charge state
			tempPosition = tempPosition - 2;
			if( sLine[ tempPosition ] == '.' )
				bIsPeptideLine = true;
		}
	}
	if( bIsPeptideLine )
		return peptideLine;
	
	/*
	 * if the fourth white-space-delimited  field is the sequence coverage field, which
	 * has a percentage sign as the last character,
	 * then this line is a protein line
	 */
	
	// extract four white-spacefields from this line
	// and put them into the the temp strings
	string sTemp1;
	string sTemp2;
	string sTemp3;
	string sTemp4;
	istringstream issStream( sLine );
	issStream >> sTemp1 >> sTemp2 >> sTemp3 >> sTemp4;


	if( sTemp4[ ( sTemp4.length() - 1 ) ] == '%' )
		return proteinLine;
	
	if( sTemp1[ ( sTemp1.length() - 2 ) ] == '.' )
		return peptideLine;

	return unknownLine;

}







