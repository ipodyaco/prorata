
#ifndef PRORATASICREADER_H
#define PRORATASICREADER_H

#include "sicReaderHelper.h"

#include <string>
#include <iostream>

#ifdef _WIN32
#include <io.h>
#define W_OK 02
#define R_OK 04
#define access( a, b ) _access( a, b )
#else
#include <unistd.h>
#endif

using namespace std;

class ProRataSICReader
{
	public:

		/* 
		 * An instance can be created either by providing the mzXML
		 * or by just creating a plain instace by default constructor
		 * and providing the filename later with setFileName function.
		 */

		ProRataSICReader( const string &sFile );
		ProRataSICReader();

		/*
		 * Perform all the clean up for parser related stuff.
		 */

		~ProRataSICReader();

		/*
		 * Run time file change is possible.
		 */

		virtual int setFileName( const string &sFile );

		/*
		 * Seem unneccessary at this moment, but might not use
		 * later.

		virtual void process();
		 */

		/*
		 * A pointer for the current identifier will be held internally
		 * so the you can call this function sequenctially to get serial
		 * access to all the Chromatograms.
		 *
		 * Note: Returns "" string if the internal pointer reaches end.
		 */

		virtual string getNextChromatogram();

		/*
		 * Reset the internal identifier pointer to the first Chromatogram
		 * identifier.
		 */

		virtual void resetIdentifierPointer()
		{	iCurrentIdentifierPointer = 0; 	}

		/*
		 * Given a Chromatogram identifier, you will get the entire 
		 * contents of the <CHROMATOGRAM> tag including it.
		 */

		virtual string getChromatogram( int iChroIdentifier ) const;

		/*
		 * Retrieve MS filename found in the first tag of SIC file.
		 * It is present as a attribute in <<EXTRACTED_ION_CHROMATOGRAMS>
		 * tag.
		 */

		virtual string getMassSpecFile() const;

		/*
		 * Retrieve the identifier of the *first* Chromatograms.
		 * It is present as a attribute in <EXTRACTED_ION_CHROMATOGRAMS>
		 * tag.
		 */

		virtual int getFirstIdentifier() const;

		/*
		 * Retrieve the identifier of the *last* Chromatograms.
		 * It is present as a attribute in <EXTRACTED_ION_CHROMATOGRAMS>
		 * tag.
		 */

		virtual int getLastIdentifier() const;

		/*
		 * Retrieve the name of the software that this file is generated 
		 * from. It is present as a attribute in <EXTRACTED_ION_CHROMATOGRAMS>
		 * tag.
		 */

		virtual string getProgramName() const;

		/*
		 * Retrieve the version of the software that this file is generated 
		 * from. It is present as a attribute in <EXTRACTED_ION_CHROMATOGRAMS>
		 * tag.
		 */

		virtual string getVersionNumber() const;

		/*
		 * Retrieve all the header information from the file. Essentially reading
		 * all the attributes of <EXTRACTED_ION_CHROMATOGRAMS> tag from the file.
		 */

		//virtual void extractHeader();

	protected:

		string sFileName;
		string sMSFile;
		
		int iFirstIdentifier;
		int iLastIdentifier;

		string sProgramName;
		string sVersionNumber;

	private:
		int iCurrentIdentifierPointer;
};

#endif //PRORATASICREADER_H

