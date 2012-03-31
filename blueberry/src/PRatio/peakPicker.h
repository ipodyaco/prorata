
#ifndef PEAKPICKER_H
#define PEAKPICKER_H

#include "chromatogram.h"
#include "smoother.h"
#include "proRataConfig.h"
#include "pca.h"

#include <vector>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <map>

using namespace std;

class PeakPicker
{
	public:
		PeakPicker();
		~PeakPicker();

		/*
		 * compute the covariance chromatogram, peak boundaries, peak shifts
		 * shift vdReferenceChroInput, if there is a peak shift detected
		 */
		bool process(   const vector< float > & vfMS2TimeInput,
				const vector< float > & vfRetentionTimeInput, 
				vector< double > & vdTreatmentChroInput,  
				vector< double > & vdReferenceChroInput
				);

		/*
		 * Various accessors to get the results.
		 */
		vector< double > & getCovarianceChro(){ return vdCovarianceChro; }
		int getLeftValley() { return iLeftValley; }
		int getRightValley() { return iRightValley; }
		bool getPeakValidity() { return bPeakValidity; }
		double getPeakHeightPPC() { return dPeakHeightPPC; }
		int getPeakShift() { return iNumberofScansShifted; }


	private:

		// find the left valley and right valley in the vdCovarianceChro
		// and then compute the PPC peak height
		void findPeak();

		// two helpers functions for findPeak()
		int findNextLeftValley( int iCurrentLeftValley );
		int findNextRightValley( int iCurrentRightValley );	

		// shift the vdReferenceChroInput by iScanShiftReference number of scans
		// compute vdCovarianceChro with the vdTreatmentChroInput and vdReferenceChroInput
		void computeCovarianceChro( int iScanShiftReference,
				vector< double > vdTreatmentChroInput,  
				vector< double > vdReferenceChroInput );

		// shift the input vdChro by the iScanShift number of scans
		// helper function for computeCovarianceChro
		void shiftChro( int iScanShift, vector< double > & vdChro );

		// Calculate the median value of the given vector
		double median( vector<double> vdData );

		// variables set by computeCovarianceChro()
		vector< double > vdCovarianceChro;
		
		// variables set by findPeak()	
		int iLeftValley;
		int iRightValley;
		double dPeakHeightPPC;
		bool bPeakValidity;
		
		// variables set by process()
		int iChroLength;
		float fLeftMS2Time;
		float fRightMS2Time;
		vector< float > vfRetentionTime; 
		int iNumberofScansShifted;
		

};

#endif //PEAKPICKER_H
