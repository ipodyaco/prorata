#include "mzReader.h"
#include <climits>

using namespace std;

mzReader::mzReader()
{
	szXMLFile = "";
	iAnalysisFirstScan = 0;
	iAnalysisLastScan = 0;
}

mzReader::~mzReader()
{
	if( szXMLFile != "" )
	{
		rampCloseFile(pFI);
	}

}

bool mzReader::setFilename( string sFilename )
{
	// Try to open the file.
	szXMLFile = sFilename;
	if ( (pFI = rampOpenFile( szXMLFile.c_str() )) == NULL)
	{
		cout << "ERROR: Could not open the given mzXML or mzData file " <<
			szXMLFile << endl;
		return false;
	}

	// Get the index offset.
	indexOffset = getIndexOffset( pFI );

	// Look in to the indexes.
	iAnalysisFirstScan = 1;
	pScanIndex = readIndex( pFI , indexOffset, &iAnalysisLastScan );

	// clear the deqMassSpecBuffer
	deqMassSpecBuffer.clear();

	return true;
}

bool mzReader::getHeaderInfo(unsigned long int iScan, 
		int * piMSLevel, double * pdPrecursorMZ, 
		int * piPeaksCount, double * pdRetentionTime)
{
	if ( iScan < iAnalysisFirstScan || iScan > iAnalysisLastScan )
	{
		return false;
	}

	// Now read the header
	readHeader( pFI, pScanIndex[iScan], &scanHeader);

	*piMSLevel = scanHeader.msLevel;
	
	// get the precursor ion's m/z for MSn scans
	if( scanHeader.msLevel > 1 )
		*pdPrecursorMZ = scanHeader.precursorMZ;
	else
		*pdPrecursorMZ = 0;
	
	*piPeaksCount = scanHeader.peaksCount;

	// convert retention time's unit from seconds to minutes
	*pdRetentionTime = ( scanHeader.retentionTime / 60 );

	return true;
}


bool mzReader::getPeaks(unsigned long int iScan, 
		int iPeaksCount, vector<float> & vfMass, 
		vector<float> & vfInten )
{
	// clear the vector
	vfMass.clear();
	vfInten.clear();
	
	// reserve the correct capacity
	vfMass.reserve( iPeaksCount );
	vfInten.reserve( iPeaksCount );
	
	int iPeaksCountCopy = iPeaksCount;
	int n = 0;
	
	if ( iScan < iAnalysisFirstScan || iScan > iAnalysisLastScan )
	{
		cout << "ERROR: Scan " << iScan << " is not between the first scan " << iAnalysisFirstScan 
			<< " and the last scan " << iAnalysisLastScan<< endl;
		return false;
	}

	if ( iPeaksCount <= 0 )
	{
	//	cout << "WARNING: Scan " << iScan << " have no peaks. " << endl;
		vfMass.push_back( 0.0 );
		vfInten.push_back( 0.0 );
		return true;
	}

	// pPeaks is an array of float dynamically allocated by RAMP
	// the float alternates between a m/z and its intensity

	float *pPeaks;

	pPeaks = readPeaks( pFI, pScanIndex[iScan] );

	float fUpperBound =  numeric_limits<float>::max() / 2.0;
	
	while ( iPeaksCountCopy-- > 0 )
	{

		/*
		vfMass.push_back( pPeaks[n] );
		n++;
		vfInten.push_back( pPeaks[n] );
		n++;
		*/
	
		// check for -1.0#QNAN
		if ( pPeaks[n] >= 0.0  && pPeaks[n] <= fUpperBound )
		{
			vfMass.push_back( pPeaks[n] );
		}
		else
		{
			cout << "WARNING: invalid m/z = " << pPeaks[n] << " at scan " << iScan << "! Re-set to zero" << endl;
			vfMass.push_back( 0.0 );
		}
		
		n++;

		if ( pPeaks[n] >= 0.0  && pPeaks[n] <= fUpperBound )
		{
			vfInten.push_back( pPeaks[n] );
		}
		else
		{
			cout << "WARNING: invalid intensity = " << pPeaks[n] << " at scan " << iScan << "! Re-set to zero" << endl;
			vfInten.push_back( 0.0 );
		}
		
		n++;
		
		
	}

	// free the memory 
	free( pPeaks );

	return true;
}

bool mzReader::getPeaksBuffered(unsigned long int iScan, 
		int iPeaksCount, vector<float> & vfMass, 
		vector<float> & vfInten )
{
	// clear the vector
	vfMass.clear();
	vfInten.clear();
	
	// check if this scan is already in the deqMassSpecBuffer
	for( iterBuffer = deqMassSpecBuffer.begin(); iterBuffer != deqMassSpecBuffer.end(); ++iterBuffer )
	{
		// if it is, get it from the buffer
		if( iterBuffer->iScan == iScan )
		{
			vfMass = iterBuffer->vfMass;
			vfInten = iterBuffer->vfInten;
			return true;
		}
	}
	
	// if this scan is not in the buffer, get it with the normal function getPeaks 
	vector< float > vfMyMass;
	vector< float > vfMyInten;
	bool bSucess = getPeaks( iScan, iPeaksCount, vfMyMass, vfMyInten );
	vfMass = vfMyMass;
	vfInten = vfMyInten;

	// remove the oldest MassSpectrum from the Buffer
	if( deqMassSpecBuffer.size() > BUFFER_SIZE )
	{
		deqMassSpecBuffer.pop_front();
	}
	// save this MassSpectrum into the Buffer
	MassSpectrum currentMassSpectrum;
	currentMassSpectrum.iScan = iScan;
	currentMassSpectrum.vfMass = vfMyMass;
	currentMassSpectrum.vfInten = vfMyInten;
	deqMassSpecBuffer.push_back( currentMassSpectrum );

	return bSucess;
}
