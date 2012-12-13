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
#include "tokenvector.h"

using namespace std;

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
		// process a peptide line and save the info into pID
		bool processPeptideLine( string sLine, Identification * pID );

		// find the number of the field that matches the targetWord
		bool findField(const TokenVector & words, string targetWord, unsigned int & iFieldNumber);

};

#endif //DTASELECTREADER_H
