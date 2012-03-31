#include "mzReader.h"
#include <climits>

using namespace std;

mzReader::mzReader()
{
	szXMLFile = "NONE";
}

mzReader::~mzReader()
{

}

bool mzReader::setFilename( string sFilename )
{
	// Try to open the file.
	szXMLFile = sFilename;
	if ( !indexFile.BuildIndex( szXMLFile ) )
	{
		cout << "ERROR: Could not open the given file " <<
			szXMLFile << endl;
		return false;
	}
	
	// clear the deqMassSpecBuffer
	deqMassSpecBuffer.clear();

	return true;
}

bool mzReader::getAllScanNumbers( vector<int> & viAllScanNumbers )
{
	viAllScanNumbers.clear();
	viAllScanNumbers = indexFile.scan_numbers();

	if ( viAllScanNumbers.size() == 0 )
	{
		cout << "no scan. empty MS file =" << szXMLFile << endl;
		return false;
	}

	return true;
}


bool mzReader::getHeaderInfo(unsigned long int iScan, 
		int * piMSLevel, double * pdPrecursorMZ, 
		int * piPeaksCount, double * pdRetentionTime)
{
	*piMSLevel = 1;
	*pdPrecursorMZ = 0;
	*piPeaksCount = 0;
	*pdRetentionTime = 0;

  vector <string> scan_data;
  istringstream input;

  if( !indexFile.GetScanText((int)iScan, scan_data) ) return false;
  
  vector <string>::iterator viter;
  int i = 0;
  
  for(viter = scan_data.begin(); viter != scan_data.end(); viter++) 
  {
	  ++i;
    TokenVector words((*viter), " \t\n\r");
    if(words.size() == 3 && words[0][0] == 'I' && words[1] == "RetentionTime") 
    { 
      istringstream input(words[2]);
      input >> *pdRetentionTime;
    }
    if(words[0][0] >= '0' && words[0][0] <= '9')
    {

	    *piPeaksCount =  scan_data.size() - i + 1;
	    break;
    }
  }

  return true;
}


bool mzReader::getPeaks(unsigned long int iScan, 
		vector<float> & vfMass, vector<float> & vfInten )
{

	// clear the vector
	vfMass.clear();
	vfInten.clear();
	
	vector <string> scan_data;
	if( !indexFile.GetScanText((int)iScan, scan_data) ) return false;

	vector <string>::iterator viter;
	float fMassTmp;
	float fIntenTmp;
	string sLine;

	for(viter = scan_data.begin(); viter != scan_data.end(); viter++) 
	{
		sLine = (*viter);
		// if a line starts with a digit
		if(sLine[0] >= '0' && sLine[0] <= '9')
		{
		    istringstream input(sLine);
		    input >> fMassTmp >> fIntenTmp;
		    vfMass.push_back(fMassTmp);
		    vfInten.push_back(fIntenTmp);
		}
	}

	return true;
}

bool mzReader::getPeaksBuffered(unsigned long int iScan, 
		vector<float> & vfMass, vector<float> & vfInten )
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
	bool bSucess = getPeaks( iScan, vfMyMass, vfMyInten );
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
