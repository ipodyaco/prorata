#include "idData.h"



IDdata::IDdata()
{
	//constructor
}

IDdata::~IDdata()
{
	// free the memory for each Identification instances
	list< Identification * >::iterator iter;
	for( iter = lpIDlist.begin(); iter != lpIDlist.end(); iter++)
	{
		delete (*iter);
	}
	
}

bool IDdata::setFilename( string sFilename )
{
	bool bSuccess;

	// pass the lpIDlist reference to DTASelectReader
	// and allow DTASelectReader to populate it
	DTASelectReader reader;
	bSuccess = reader.getIDlist( sFilename, lpIDlist );

	// set the iFirstMS2 and iLastMS2 in each Identification
	list< Identification * >::iterator iter;
	for( iter = lpIDlist.begin(); iter != lpIDlist.end(); ++iter)
	{
		setFirstAndLastMS2( *iter );
	}
	
	return bSuccess;
}


	

	
void IDdata::consolidateIDlist( MSdata *pMSdata )
{
	string sBaseFilename = pMSdata->getBaseFilename();

	/*
	 * compare every possible pairs of Identifications in lpIDlist, (*iter1) and (*iter0)
	 * if they can be merged, merge (*iter1) into (*iter0) and then mark (*iter0) as to be deleted
	 * by changing its charge state to -9999.
	 */
	list< Identification * >::iterator iter0;
	list< Identification * >::iterator iter1;

	setRetentionTime( pMSdata );
	

	// a vector containning the iterators for the Identifications that have to be deleted
	vector< list< Identification * >::iterator > vIterDeletion;
	
	for( iter0 = lpIDlist.begin(); iter0 != lpIDlist.end(); ++iter0 )
	{
		// determine if the current (*iter0) is marked as to be deleted
		// if so, skip it
		if( (*iter0)->iChargeState == -9999 )
			continue;

		// determine if this ID is from this MS file
		if( !isFilenameMatched( sBaseFilename, (*iter0) ) )
			continue;

		// initialize iter1 as the one immediately after iter0
		iter1 = iter0;
		for( ++iter1; iter1 != lpIDlist.end(); ++iter1 )
		{
			// determine if the current (*iter1) is marked as to be deleted
			// if so, skip it
			if( (*iter1)->iChargeState == -9999 )
				continue;

			// determine if they are of the same charge state
			// if not, skip it
			if( (*iter0)->iChargeState != (*iter1)->iChargeState )
				continue;

			// determine if they are of the same sequence
			// if not, skip it
			if( (*iter0)->sSequence != (*iter1)->sSequence )
				continue;			

			// determine if this ID is from this MS file
			if( !isFilenameMatched( sBaseFilename, (*iter1) ) )
				continue;

			// determine if they have adajent MS2 scans
			if( isMS2Adjacent( (*iter0), (*iter1) ) )
			{
				// if so, merge (*iter1)'s info into (*iter0)
				mergeID( (*iter0), (*iter1) );

				// mark (*iter1) as to be deleted by changing its charge state
				// to -9999
				(*iter1)->iChargeState = -9999;

				// save iter1 to the vector of iterators to be deleted
				vIterDeletion.push_back( iter1 );
			}
		}
	}

	// for all iterators that need to be deleted
	// first free the memory and then erase the iterator from lpIDlist
	// Deletions in list only invalidate iterators and references to the deleted elements. 
	for(unsigned int i = 0; i< vIterDeletion.size(); ++i )
	{
		iter1 = vIterDeletion[i];
		delete (*iter1);
		lpIDlist.erase(iter1);
	}

}



bool IDdata::getIDvector( string sBaseFilename, vector< Identification * > & vpIDvector )
{

	// clear the vector
	vpIDvector.clear();
	
	list< Identification * >::iterator iter;
	
	// save to the vpIDvector all the Identifications matching
	// the given Base Filename
	for( iter = lpIDlist.begin(); iter != lpIDlist.end(); ++iter)
	{
		if( isFilenameMatched( sBaseFilename, (*iter) ) )
			vpIDvector.push_back( (*iter) );

	}

	// creat the function object
	lessID lessIDsort;

	// use the genetic algorithm, "sort", to sort vpIDvector
	sort( vpIDvector.begin(), vpIDvector.end(), lessIDsort );

	if( vpIDvector.size() < 1 )
		return false;
	else
		return true;

}

void IDdata::showMeAll()
{
	// output everything to a temp file
	ofstream logFile( "logTemp.txt" );
	Identification * pID1;
	vector< MS2Scoring > vMS2Scoring;
	vector< Protein > vProtein;
	list< Identification * >::iterator iter;
	unsigned int i;
	logFile << " Totol Identification Number = " << lpIDlist.size () << endl;

	// print out everything in lpIDlist
	for( iter = lpIDlist.begin(); iter != lpIDlist.end(); iter++)
	{
		pID1 = (*iter);
		logFile << "Charge state =  " << pID1->iChargeState << " first MS2 = " << pID1->vMS2Scoring[pID1->iFirstMS2].iMSMSscan << " last MS2 = " << pID1->vMS2Scoring[pID1->iLastMS2].iMSMSscan <<  "	Sequence =  " << pID1->sSequence  << endl;
		vMS2Scoring = pID1->vMS2Scoring;
		vProtein = pID1->vProtein;
		for( i = 0; i < vMS2Scoring.size(); ++i)
			logFile << "	Score =  " << vMS2Scoring[i].fScore << "	iMSMSscan =  " << vMS2Scoring[i].iMSMSscan  << "	fRT = " << vMS2Scoring[i].fRetentionTime << "	sIDfilename =	" <<  vMS2Scoring[i].sIDfilename << endl;
		for( i = 0; i < vProtein.size(); ++i )
			logFile << "		Locus =  " << vProtein[i].sLocus << "	Description =  " << vProtein[i].sDescription << endl;
		
	}

}

bool IDdata::isFilenameMatched( string sBaseFilename, const Identification * pID )
{
	// if the base filename is a part of all MS2Scoring's sIDfilename 
	// in the vMS2Scoring vector, then the fileaname is matched
	int iMS2Count = pID->vMS2Scoring.size();
	string sIDfilename;
	for( int i = 0; i < iMS2Count; ++i )
	{
		sIDfilename = pID->vMS2Scoring[i].sIDfilename;
		if( sIDfilename.find( sBaseFilename, 0 ) == string::npos )
		{
			return false;
		}
	}
	return true;
}

bool IDdata::isMS2Adjacent( Identification * pID0, Identification * pID1 )
{
	// determine if these two IDs have adjacent MS2 scans
	// the RT range for ID0 is between fLowerRT0 and fUpperRT0
	// then expand it by the allowing retention time window between them
	float fLowerRT0 = pID0->vMS2Scoring[pID0->iFirstMS2].fRetentionTime - ProRataConfig::getMinutesBetweenMS2();
	float fUpperRT0 = pID0->vMS2Scoring[pID0->iLastMS2].fRetentionTime + ProRataConfig::getMinutesBetweenMS2();

	// the RT range for ID1 is between  fLowerRT1 and fUpperRT1
	float fLowerRT1 = pID1->vMS2Scoring[pID1->iFirstMS2].fRetentionTime;
	float fUpperRT1 = pID1->vMS2Scoring[pID1->iLastMS2].fRetentionTime;
	
	// determine if the two RT ranges overlap
	if( ( fLowerRT0 <= fUpperRT1 ) && ( fLowerRT1 <= fUpperRT0 ) )
		return true;
	else
		return false;
}

void IDdata::mergeID( Identification * pID0, Identification * pID1 )
{
	unsigned int i0;
	unsigned int i1;
	bool isMatched = false;

	// merge Protein vector
	// if ID1 has a Protein that is not in ID0's vProtein
	// then save that Protein to ID0's vProtein 
	for( i1 = 0; i1 < pID1->vProtein.size(); ++i1)
	{
		isMatched = false;
		for( i0 = 0; i0 < pID0->vProtein.size(); ++i0)
		{
			// if the two Protein intances have the same locus
			// then they are the same
			if( pID0->vProtein[i0].sLocus == pID1->vProtein[i1].sLocus )
			{
				isMatched = true;
				break;
			}
		}
		if( !isMatched )
			pID0->vProtein.push_back( pID1->vProtein[i1] );
	}	

	// merge MS2Scoring vector
	// if ID1 has a MS2Scoring are not in ID0's vMS2Scoring
	// then save it to  ID0's vMS2Scoring
	for( i1 = 0; i1 < pID1->vMS2Scoring.size(); ++i1)
	{
		isMatched = false;
		for( i0 = 0; i0 < pID0->vMS2Scoring.size(); ++i0)
		{
			// if the two MS2Scoring have the same MS2 scan number
			// and the same ID filename, then they are the same
			if( pID0->vMS2Scoring[i0].iMSMSscan == pID1->vMS2Scoring[i1].iMSMSscan &&
					pID0->vMS2Scoring[i0].sIDfilename == pID1->vMS2Scoring[i1].sIDfilename )
			{
				isMatched = true;
				break;
			}
		}
		if( !isMatched )
			pID0->vMS2Scoring.push_back( pID1->vMS2Scoring[i1] );
	}

	// re-set the iFirstMS2 and iLastMS2 for ID0
	setFirstAndLastMS2(pID0);
}


void IDdata::setFirstAndLastMS2( Identification * pID )
{
	// abort, if this Identification has zero MS2Scoring
	if( pID->vMS2Scoring.size() < 1 )
		return;
	
	// initialize iFirstMS2ScanNumber and iLastMS2ScanNumber to
	// the first MS2Scoring's MS2 scan number
	int iFirstMS2 = 0;
	int iLastMS2 = 0;
	unsigned long int iFirstMS2ScanNumber = pID->vMS2Scoring[iFirstMS2].iMSMSscan;
	unsigned long int iLastMS2ScanNumber = pID->vMS2Scoring[iLastMS2].iMSMSscan;
	for(unsigned int i = 1; i < pID->vMS2Scoring.size(); ++i )
	{
		if( iFirstMS2ScanNumber > pID->vMS2Scoring[i].iMSMSscan )
		{
			iFirstMS2 = i;
			iFirstMS2ScanNumber = pID->vMS2Scoring[iFirstMS2].iMSMSscan;
		}
		
		if( iLastMS2ScanNumber < pID->vMS2Scoring[i].iMSMSscan )
		{
			iLastMS2 = i;
			iLastMS2ScanNumber = pID->vMS2Scoring[iLastMS2].iMSMSscan;
		}
	}

	// save results back to pID
	pID->iFirstMS2 = iFirstMS2;
	pID->iLastMS2 = iLastMS2;
	
}


void IDdata::setRetentionTime( MSdata *pMSdata )
{
	list< Identification * >::iterator iter;
	
	unsigned int i;
	unsigned long int iScanNumber;
	float fRetentionTime;
	for( iter = lpIDlist.begin(); iter != lpIDlist.end(); iter++)
	{
		if( isFilenameMatched( pMSdata->getBaseFilename(), (*iter) ) )
		{
			for( i = 0; i < (*iter)->vMS2Scoring.size(); i++ )
			{
				iScanNumber = (*iter)->vMS2Scoring[i].iMSMSscan;
				if( pMSdata->getTime4Scan( iScanNumber, fRetentionTime ) )
					(*iter)->vMS2Scoring[i].fRetentionTime =  fRetentionTime;
			}
		}

	}
}



