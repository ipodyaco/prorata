#ifndef DTASELECTREADER_H
#define DTASELECTREADER_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream> 
#include <vector>
#include <list>
#include "chromatogram.h"
#include "idData.h"

#define CSTRING_SIZE 256


using namespace std;

// enumerate the line types in a DTASelect-filter file
enum LineType
{
	// a peptide line
	peptideLine,

	// a protein line
	proteinLine,

	// the line proceeding the peptide/protein lines
	startLine,

	// the line trailing the peptide/protein lines
	endLine,

	// the unknown line
	unknownLine
};

class DTASelectReader
{
	public:
		DTASelectReader();
		~DTASelectReader();

		/* 
		 * two overloaded functions for get a ID list from a DTASelect-filter file
		 * the ID list is returned either by value or by reference
		 */
		list< Identification * > getIDlist( string sFilename );
		bool getIDlist( string sFilename, list< Identification* > & lpIDlist );

	private:

		// process a protein line and save the info into vProtein
		bool processProteinLine( string sLine, vector< Protein > & vProtein );

		// process a peptide line and save the info into pID
		bool processPeptideLine( string sLine, const vector< Protein > & vProtein, Identification * pID );

		// determine the line type of the current line
		LineType whatIsThisLine( string sLine );
		

};

#endif //DTASELECTREADER_H
