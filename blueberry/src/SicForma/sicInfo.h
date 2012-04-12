#ifndef SICINFO_H
#define SICINFO_H

#include <string>
#include <vector>
#include <sys/stat.h> 
#include <sys/types.h> 
#include "chromatogram.h"
#include "directoryStructure.h"
#include "idData.h"
#include "msData.h"
#include "isotopologue.h"
#include "peptideInfo.h"
#include "peptideRatio.h"
#include "peptideLabelFree.h"


using namespace std;
// this is the number of the chromatogram in a xic file
// const int CHRO_COUNT = 100;
class SICinfo
{
	public:
		SICinfo();
		~SICinfo();

		/*
		 * set the filenames for the ID files and the MS files, which
		 * is stored in sIDfilename and vsMZfilename
		 * the ID filename is sIDfilenameInput
		 * the MS filenames are all files with extension name of mzXML/mzData
		 * in the working directory.
		 */
		bool setFilename( string sIDfilenameInput );

		/*
		 * instantiate the IDdata, MSdata and Isotopologue classes
		 * creat extract chromatograms and write them into the file one by one
		 */
		bool process(vector< PeptideInfo * > & vpPeptideInfo);
		
		
	private:

		/*
		 * given an Identification, a chromatogram is extracted from the current MS file that pMSdata is pointing to 
		 * its identifier is the parameter iIdentifier, the chromatogram is returned by reference
		 */
		bool extractChromatogram( const Identification & idInput, string sMSfilenameInput, int iIdentifier, Chromatogram & chroOutput );

		/*
		 * compose a Root element start tag for the output file, sXicFilename
		 * and then write it to the file.
		 */
		bool writeRootElementStartTag( string sXicFilename, int iFirstIdentifier, int iLastIdentifier );
		
		/* 
		 * ID file's filename
		 */
		string sIDfilename;

		/*
		 * the pointer to the ID data object for processing the given ID file
		 */
		IDdata* pIDdata;

		/*
		 * the mzXML/mzData filenames in the working directory
		 */
		vector< string > vsMZfilename;

		/*
		 * a pointer to an MS data object
		 * during process(), this pointer is moved down the list of MZ files
		 */
		MSdata * pMSdata;

		/* 
		 * a pointers to Isotopologues objects
		 */
		vector< Isotopologue * > vpIsotopologue;

};

#endif //SICINFO_H 
