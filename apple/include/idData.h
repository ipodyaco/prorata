#ifndef IDDATA_H
#define IDDATA_H

#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>
#include "chromatogram.h"
#include "dtaSelectReader.h"
#include "msData.h"

using namespace std;

// forward declaration
class DTASelectReader;
/*
 * a function object, which is passed to the generic algorithm sort() to sort
 * a vector of Identification pointers. the Identification pointers are sorted by
 * their first MSMS scan number
 */
class lessID
{
	public:
		// overload the operator ()
		// return true if the MSMS scan number of ID1 is less then that of ID2.
		// return false, otherwise.
		bool operator() ( const Identification * pID1, const Identification * pID2 ) const
		{
			unsigned long int iMSMSscanID1 = pID1->vMS2Scoring[pID1->iFirstMS2].iMSMSscan;
			unsigned long int iMSMSscanID2 = pID2->vMS2Scoring[pID2->iFirstMS2].iMSMSscan;
			if( iMSMSscanID1 < iMSMSscanID2 )
				return true;
			else
				return false;
		}
};

/*
 * a class for holding all Identification data in a list
 * it should be used by calling functions in the following steps
 * 1) setFilename():  populate the list with DTASelectReader or other ID readers
 * 2) setRetentionTime():  set the retention time for each MSMS scan
 * 3) consolidateIDlist():  consolidate the list by removing redundance identification and merging adjacent MS2 scans
 * 4) getIDvector(): retrieve a subset of the list in the sorted order
 */
class IDdata
{
	public:
		IDdata();
		
		/*
		 * free memory
		 */
		~IDdata();
		
		/*
		 * set the filename for a MS/MS identification results
		 * and select appropriet ID reader to populate the Identification list
		 */
		bool setFilename( string sFilename );

		/*
		 * consolidate the identification list by merging identifications for
		 * the same peptide at the same charge state identified at about the same time from the same MS file
		 * the merged Identifications most likely come from the different isotopologues of the same peptide 
		 */
		void consolidateIDlist( MSdata * pMSdata );

		/*
		 * get a vector of identifications from a mzXML/mzData file
		 * the identification is sorted by their first MS/MS scan number
		 */
		bool getIDvector( string sBaseFilename, vector< Identification * > & vpIDvector );

		/*
		 * debug function
		 * print out all Identification to a temp file logTemp.txt
		 */
		void showMeAll();
		
	private:
		// a list containing all ID instances's pointer
		list< Identification * > lpIDlist;
		
		// set the iFirstMS2 and iLastMS2 for an Identification
		void setFirstAndLastMS2( Identification * pID );
		
		// is this ID's filename matching the BaseFile name?
		bool isFilenameMatched( string sBaseFilename, const Identification * pID );

		// are the two Identifications from the same mzXML file and close in their retention time?
		bool isMS2Adjacent( Identification * pID0, Identification * pID1);

		// merge the two IDs
		void mergeID( Identification * pID0, Identification * pID1 );
		
		void setRetentionTime( MSdata *pMSdata );
	
};


#endif //IDDATA_H
