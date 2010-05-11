
#ifndef PRORATAMASSSPECTRUMDATA_H
#define PRORATAMASSSPECTRUMDATA_H

#include "peptideRatio.h"
#include "directoryStructure.h"
#include "mzReader.h"

class ProRataMassSpectrumData
{
	public:
		ProRataMassSpectrumData();
		~ProRataMassSpectrumData();

		bool getScan( string sFilename, unsigned long int iScan, 
				vector< double > & vdMZ, vector< double > & vdRelativeIntensity );
		bool getScan( string sFilename, unsigned long int iScan, 
				vector< double > & vdMZ, vector< double > & vdRelativeIntensity, double & dPrecursorMZ );
		
	private:
		map<string, mzReader*> mapReaderFileMappings;

};
#endif //PRORATAMASSSPECTRUMDATA_H
