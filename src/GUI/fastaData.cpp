#include "fastaData.h"

FASTAdata::FASTAdata()
{
	// constructor
}

FASTAdata::FASTAdata( string sFilename )
{
	readFASTA( sFilename );
}

FASTAdata::~FASTAdata()
{
	// destructor
}

bool FASTAdata::readFASTA( string sFilename )
{
	ifstream fileFASTA( sFilename.c_str() );
	string sCurrentLine( "" );

	string sTempSequence( "" );
	string sCurrentLocus( "" );
	
	if( !fileFASTA )
	{
		cout << "Cannot open FASTA file: " << sFilename << endl;
		return false;
	}

	// skip to the first description line
	while( getline(fileFASTA, sCurrentLine) )
	{
		if( isDesciptionLine( sCurrentLine ) )
			break;
	}

	// process the lines
	do{
		if( isDesciptionLine( sCurrentLine ) )
		{
			saveSequence( sCurrentLocus, sTempSequence );
			sTempSequence = "";
			sCurrentLocus = processDescriptionLine( sCurrentLine );
		}
		else
		{
			sTempSequence = processSequenceLine( sTempSequence, sCurrentLine );
		}
	}
	while( getline(fileFASTA, sCurrentLine) );

	return true;

}

string FASTAdata::getDescription( string sLocus )
{
	string sDescription;
	map< string, string >::iterator iter = mLocus2Description.find( sLocus );
	
	if( iter == mLocus2Description.end() )
	{
		sDescription = "";
	}
	else
	{
		sDescription = (*iter).second;
	}
	
	return sDescription;
}

string FASTAdata::getProteinSequence( string sLocus )
{
	string sProteinSequence;
	map< string, string >::iterator iter = mLocus2Sequence.find( sLocus );
	
	if( iter == mLocus2Sequence.end() )
	{
		sProteinSequence = "";
	}
	else
	{
		sProteinSequence = (*iter).second;
	}
	
	return sProteinSequence;
}

void FASTAdata::computeDistances( string sLocus, string sPeptideSequence, int *piNterminalDistance, int *piCterminalDistance)
{
	string sProteinSequence;
	map< string, string >::iterator iter = mLocus2Sequence.find( sLocus );
	
	// exit if the locus cannot be found or the input peptide string is too short
	if( iter == mLocus2Sequence.end() || sPeptideSequence.length() < 3 )
	{
		*piNterminalDistance = 0;
		*piCterminalDistance = 0;
		return;
	}
	
	sProteinSequence = (*iter).second;

	// extract sequence from X.SEQ.X format
	int iEndingPeriod = (sPeptideSequence.length() - 2);
	if( sPeptideSequence[1] == '.' && sPeptideSequence[iEndingPeriod] == '.' )
	{
		sPeptideSequence = sPeptideSequence.substr( 2, (iEndingPeriod - 2) );
	}

	// format the peptide sequence
	sPeptideSequence = formatSequence( sPeptideSequence );

	if( sProteinSequence.find( sPeptideSequence ) != string::npos )
	{
		// adding 1 is because the protein sequence is indexed from 1
		*piNterminalDistance = sProteinSequence.find( sPeptideSequence ) + 1 ;
		*piCterminalDistance = sProteinSequence.length() - sPeptideSequence.length() - *piNterminalDistance + 2 ;	
	}
	else
	{
		*piNterminalDistance = 0;
		*piCterminalDistance = 0;
	}
	
}

void FASTAdata::computeCoveragePlot( const string & sProteinSequence, string sPeptideSequence, int & iNtermCord, int & iCtermCord )
{
	sPeptideSequence = formatSequence( sPeptideSequence );
	string::size_type pos = sProteinSequence.find( sPeptideSequence );
	if( pos != string::npos )
	{
		// adding 1 is because the protein sequence is indexed from 1
		iNtermCord = pos + 1 ;
		iCtermCord = iNtermCord + sPeptideSequence.length();	
	}
	else
	{
		iNtermCord = 0;
		iCtermCord = 0;
	}
}


bool FASTAdata::isDesciptionLine( string sLine )
{
	char cFirstCharacter = sLine[0];
	if( cFirstCharacter == '>' )
		return true;
	else
		return false;
}

string FASTAdata::processDescriptionLine(string sDescriptionLine)
{
	string sLocus;
	string sSecondWord;
	string sDescription;
	istringstream inputStream( sDescriptionLine );
	
	// get the first word as the locus
	// get the second word as the begin of the description
	inputStream >> sLocus >> sSecondWord;
	
	// remove the '>'
	sLocus.erase(0, 1);
	
	// get the rest of string as the description
	int iBeginDescription = sDescriptionLine.find( sSecondWord );
	sDescription = sDescriptionLine.erase( 0, iBeginDescription );

	// save the description
	mLocus2Description[ sLocus ] = sDescription;

	return sLocus;
	
}

string FASTAdata::processSequenceLine( string sTempSequence, string sSequenceLine )
{
	// determine if sSequenceLine is a comment line, which should begin with ';'
	// if it is, do nothing to sTempSequence
	// if it isn't, append sSequenceLine to sTempSequence
	char cFirstCharacter = sSequenceLine[0];
	if( cFirstCharacter == ';' )
		return sTempSequence;
	else
		return sTempSequence.append( sSequenceLine );
}

void FASTAdata::saveSequence( string sLocus, string sSequence )
{
	mLocus2Sequence[ sLocus ] = formatSequence( sSequence );
}

string FASTAdata::formatSequence( string sSequence )
{
	int iLength = sSequence.length();
	// check if this sequence is of the format X.XXXX.X
	// its minimum length is 5
	if(iLength > 4 )
	{
		// its second residues from both the left side and the right side should be '.'
		if( sSequence[1] == '.' && sSequence[ ( iLength - 2 ) ] == '.')
		{
			// if so, extract the sequence between '.'
			sSequence = sSequence.substr( 2, (iLength - 4) );
		}
	}
	
	unsigned int i = 0;
	while( i < sSequence.length() )
	{
		if( isalpha( sSequence[i] ) )
		{
			// change the lower case letters to the upper case letters
			sSequence[i]=toupper(sSequence[i]);
			i++;
		}
		else
		{
			// remove the char that is not a letter of the alphabet
			// and i stay where it is
			sSequence.erase( i, 1 );
		}
	}
	return sSequence;
}
