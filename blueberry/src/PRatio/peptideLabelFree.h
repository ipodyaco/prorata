
#ifndef PEPTIDELABELFREE_H
#define PEPTIDELABELFREE_H

//#define DEBUG

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <math.h>

#include "chromatogram.h"
#include "tinyxml.h"
#include "peakPicker.h"
#include "smoother.h"
#include "peptideInfo.h"
#include "proRataConfig.h"

using namespace std;

class PeptideLabelFree
{
	public:
		PeptideLabelFree();
		~PeptideLabelFree();

		/*
		 * read a chromatogram from the file
		 * detect the chromatographic peak
		 * compute the abundance ratio
		 * compute the signal-to-noise ratio
		 */
		bool process( const Chromatogram & inputChro );

		// accessors
		int getIdentifier(){ return peptideChro.getIdentifier(); }
		string getSequence(){ return peptideChro.getSequence(); }
		int getChargeState(){ return peptideChro.getChargeState(); }
		float getMaximumScore(){ return peptideChro.getMaximumScore(); }
		bool getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription );
		
		double getPeakHeight() { return dPeakHeight; }
		double getPeakArea() { return dPeakArea; }
		double getPeakSNR() { return dReferencePeakSNR; }
		
		bool getValidity(){ return bValidity; }
		
		// the accessor functions for ploting mass spec in GUI
		string getMSfilename(){ return peptideChro.getMSfilename(); }
		bool getReferenceMZrange( vector< float > & vfLowerMZ, vector< float > & vfUpperMZ );

		// data for drawing selected ion chromatograms in GUI
		const vector< double > & getReferenceChro(){ return vdReferenceChro; }
		const vector< float > & getRetentionTime(){ return vfRetentionTime; }
		const vector< float > & getMS2Time(){ return vfMS2Time; }
		vector< unsigned long int > getMS2ScanNumber(){ return peptideChro.getMS2ScanNumber(); }
		float getLeftValleyTime(){ return vfRetentionTime[ iLeftValley ]; }
		float getRightValleyTime(){ return vfRetentionTime[ iRightValley ]; }

		// data for drawing PCA plot in GUI
		const vector< double > & getReferencePeak(){ return vdReferencePeak; }
		vector< unsigned long int > getScanNumberPeak();
		
		// data for text area
		vector< string > getAllIDfilename(){ return peptideChro.getAllIDfilename(); }
		unsigned long int getFirstScanNumber(){ return peptideChro.getScanVector().front(); }
		unsigned long int getLastScanNumber(){ return peptideChro.getScanVector().back(); }
		
	private:

		bool computeRatio();

		double logBase2( double dNum )
		{ return (double) ( log( (double)dNum ) / log( (double)2 ) ); }
		
		bool checkChromatogramValidity();
		
		double computePeakSNR( vector< double > vdChro );

		bool quantifyWithPeak();
		
		double dPeakHeight;
		double dPeakArea;
		double dReferencePeakSNR;
		bool bValidity;
		
		PeakPicker peptidePeakPicker;
		Chromatogram peptideChro;
		
		string sChroName;
		
		vector< double > vdReferenceChro;
		vector< double > vdReferencePeak;

		int iLeftValley;
		int iRightValley;
		
		vector< float > vfRetentionTime;
		vector< unsigned long int > viScanVector;
		vector< float > vfMS2Time;
};


#endif //PEPTIDELABELFREE_H
