
#ifndef PEPTIDERATIO_H
#define PEPTIDERATIO_H

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
#include "pca.h"
#include "smoother.h"
#include "peptideInfo.h"
#include "proRataConfig.h"

using namespace std;

class PeptideRatio
{
	public:
		PeptideRatio();
		~PeptideRatio();

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
		
		double getPeakHeightRatio() { return dPeakHeightRatio; }
		double getPeakAreaRatio() { return dPeakAreaRatio; }
		double getTreatmentPeakSNR() { return dTreatmentPeakSNR; }
		double getReferencePeakSNR() { return dReferencePeakSNR; }
		double getCovariancePeakSNR() { return computePeakSNR(peptidePeakPicker.getCovarianceChro()); }
		double getPCARatio() { return dPCARatio; }
		double getPCASN() { return dPCASN; }
		void getEngenvectors( double & dPC1x, double & dPC1y, double & dPC2x, double & dPC2y );
		void getEngenvalues( double & dPC1EV, double & dPC2EV ); 
		
		double getLog2PeakHeightRatio() { return logBase2(dPeakHeightRatio) ; }
		double getLog2PeakAreaRatio() { return logBase2(dPeakAreaRatio); }
		double getLog2PCARatio() { return logBase2(dPCARatio); }
		double getLog2PCASN() { return logBase2(dPCASN); }
		
		bool getValidity(){ return bValidity; }
		
		int getPeakShift(){ return peptidePeakPicker.getPeakShift(); }

		// the accessor functions for ploting mass spec in GUI
		unsigned long int getFullScan4Time( float fTime ){ return peptideChro.getFullScan4Time( fTime ); }
		unsigned long int getNextFullScan( unsigned long int iCurrentFullScan, bool bMoveForward = true );
	//	float getFullScanTime4Time( float fTime ){ return peptideChro.getFullScanTime4Time( fTime ); }
		float getFullScanTime( unsigned long int iScan);
		string getMSfilename(){ return peptideChro.getMSfilename(); }
		bool getReferenceMZrange( vector< float > & vfLowerMZ, vector< float > & vfUpperMZ );
		bool getTreatmentMZrange( vector< float > & vfLowerMZ, vector< float > & vfUpperMZ );

		// data for drawing selected ion chromatograms in GUI
		const vector< double > & getReferenceChro(){ return vdReferenceChro; }
		const vector< double > & getTreatmentChro(){ return vdTreatmentChro; }
		const vector< double > getCovarianceChro(){ return peptidePeakPicker.getCovarianceChro(); }
		const vector< float > & getRetentionTime(){ return vfRetentionTime; }
		const vector< float > & getMS2Time(){ return vfMS2Time; }
		vector< unsigned long int > getMS2ScanNumber(){ return peptideChro.getMS2ScanNumber(); }
		float getLeftValleyTime(){ return vfRetentionTime[ iLeftValley ]; }
		float getRightValleyTime(){ return vfRetentionTime[ iRightValley ]; }

		// data for drawing PCA plot in GUI
		const vector< double > & getReferencePeak(){ return vdReferencePeak; }
		const vector< double > & getTreatmentPeak(){ return vdTreatmentPeak; }
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
		
		bool quantifyWithPCA();

		double computePeakSNR( vector< double > vdChro );

		bool quantifyWithPeak();
		
		double dPeakHeightRatio;
		double dPeakAreaRatio;
		double dReferencePeakSNR;
		double dTreatmentPeakSNR;
		double dPCARatio;
		double dPCASN;
		bool bValidity;
		
		double dxCord1; 
		double dyCord1;	
		double dxCord2;
		double dyCord2;

		double dFirstEigenValue;
		double dSecondEigenValue;
		
		PeakPicker peptidePeakPicker;
		Chromatogram peptideChro;
		
		vector< double > vdReferenceChro;
		vector< double > vdTreatmentChro;
		
		vector< double > vdReferencePeak;
		vector< double > vdTreatmentPeak;


		int iLeftValley;
		int iRightValley;
		
		vector< float > vfRetentionTime;
		vector< unsigned long int > viScanVector;
		vector< float > vfMS2Time;
};


#endif //PEPTIDERATIO_H
