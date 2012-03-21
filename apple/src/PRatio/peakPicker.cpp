#include "peakPicker.h"

PeakPicker::PeakPicker()
{
	iRightValley = 0;
	iLeftValley = 0;
	dPeakHeightPPC = 0;
	bPeakValidity = false;
	iChroLength = 0;
	fLeftMS2Time = 0;
	fRightMS2Time = 0;
	iNumberofScansShifted = 0;
}

PeakPicker::~PeakPicker()
{
	// destructor
}

bool PeakPicker::process(   const vector< float > & vfMS2TimeInput,
		const vector< float > & vfRetentionTimeInput, 
		vector< double > & vdTreatmentChroInput,  
		vector< double > & vdReferenceChroInput
		)
{
	iChroLength = vfRetentionTimeInput.size();
	
	if ( iChroLength <= ProRataConfig::getPeakPickerWindowSize() )
	{
		cout << "ERROR! The input chromatograms have too few MS1 scans!" << endl;
		return false;
	}

	if( vdTreatmentChroInput.size() != iChroLength || vdReferenceChroInput.size() != iChroLength )
	{
		cout << "ERROR! the input chromatograms are in different length!" << endl;
		return false;
	}

	if ( ( ProRataConfig::getPeakPickerWindowSize() % 2) == 0 )
	{
		cout << "ERROR! The window size for Sav-Gol smoothing has to be odd! " << endl; 
		return false;
	}
	
	if (  ProRataConfig::getPeakPickerWindowSize() < 3 )
	{
		cout << "ERROR! The window size for Sav-Gol smoothing is too small! " << endl; 
		return false;
	}
	
	vfRetentionTime = vfRetentionTimeInput; 

	// set the earliest MS2 time and the last MS2 scan time
	// used for initialize left valley and right valley
	fLeftMS2Time = *( min_element( vfMS2TimeInput.begin(), vfMS2TimeInput.end() ) ); 
	fRightMS2Time = *( max_element( vfMS2TimeInput.begin(), vfMS2TimeInput.end() ) ); 

	// the fLeftMS2Time and fRightMS2Time have to be within the RT range
	if( fLeftMS2Time < vfRetentionTime.front() )
		fLeftMS2Time = vfRetentionTime.front();

	if( fRightMS2Time > vfRetentionTime.back() )
		fRightMS2Time = vfRetentionTime.back();
	

	if( (ProRataConfig::getLeftPeakShift() < 0.001) && (ProRataConfig::getRightPeakShift() < 0.001) )
	{
		iNumberofScansShifted = 0;
		// compute the final vdCovarianceChro
		computeCovarianceChro( iNumberofScansShifted,  vdTreatmentChroInput, vdReferenceChroInput );
		// compute the final peak
		findPeak();
		return true;
	}

	/*
	 *  detect peak shift
	 */
	
	int i;

	// the calcuated peak height in PPC for all tested peak shift
	// vector< double > vdPeakHeight;
	// vector< int > viScanShift;
	map< int, double > mScanShift2Height;

	// scans that shifts left is negative and scans that shifts right is positive
	for( i = 0; i <  vfRetentionTime.size() - 3; i++ )
	{
		if( (vfRetentionTime[i] - vfRetentionTime[0]) >= ProRataConfig::getLeftPeakShift() - 0.001)
			break;
	}
	int iScanShiftLeft = -i;

	for( i = (vfRetentionTime.size() - 1); i > 2; i-- )
	{
		if( (vfRetentionTime.back() - vfRetentionTime[i]) >= ProRataConfig::getRightPeakShift() - 0.001)
			break;
	}

	int iScanShiftRight = vfRetentionTime.size() - i - 1;

	double dMaxPeakHeight = -100;
	for( i = iScanShiftLeft; i < ( iScanShiftRight + 1 ); ++i )
	{
		// viScanShift.push_back( i );
		computeCovarianceChro( i,  vdTreatmentChroInput, vdReferenceChroInput );
		findPeak();
		mScanShift2Height[i] = dPeakHeightPPC;
		if( dPeakHeightPPC > dMaxPeakHeight )
			dMaxPeakHeight = dPeakHeightPPC;
		// vdPeakHeight.push_back( dPeakHeightPPC );
	}

	if( mScanShift2Height[ 0 ] > dMaxPeakHeight*0.9 )
	{
		iNumberofScansShifted = 0;
	}
	else
	{
		iNumberofScansShifted = 0;
		bool bIsMaxFound = false;
		if( abs( iScanShiftLeft ) > iScanShiftRight )
		{
			for( i = 0; i > ( iScanShiftLeft - 1 ); --i )
			{
				if( mScanShift2Height[i] > dMaxPeakHeight*0.95 )
				{
					iNumberofScansShifted = i;
					bIsMaxFound = true;
					break;
				}
				
			}
			if( !bIsMaxFound )
			{
				for( i = 0; i < iScanShiftRight + 1 ; ++i )
				{
					if( mScanShift2Height[i] > dMaxPeakHeight*0.95 )
					{
						iNumberofScansShifted = i;
						break;
					}
				}

			}


		}
		else
		{
			for( i = 0; i < iScanShiftRight + 1 ; ++i )
			{
				if( mScanShift2Height[i] > dMaxPeakHeight*0.95 )
				{
					iNumberofScansShifted = i;
					bIsMaxFound = true;
					break;
				}
				
			}
			if( !bIsMaxFound )
			{
				for( i = 0; i > ( iScanShiftLeft - 1 ); --i )
				{
					if( mScanShift2Height[i] > dMaxPeakHeight*0.95 )
					{
						iNumberofScansShifted = i;
						break;
					}
				}

			}

		}
	}


	
	// compute the final vdCovarianceChro
	computeCovarianceChro( iNumberofScansShifted,  vdTreatmentChroInput, vdReferenceChroInput );

	// compute the final peak
	findPeak();


	
	// actually change the input vdReferenceChroInput
	shiftChro( iNumberofScansShifted, vdReferenceChroInput );

	// cout << "bPeakValidity = " <<  boolalpha <<  bPeakValidity << endl;


	
	return true;
}

void PeakPicker::computeCovarianceChro( int iScanShiftReference,
		vector< double > vdTreatmentChroInput,  
		vector< double > vdReferenceChroInput )
{
	vdCovarianceChro.clear();
	vdCovarianceChro.resize( iChroLength, 0.0 );
	
	// ensure all input chromatograms have the right length
	if( vdTreatmentChroInput.size() != iChroLength || vdReferenceChroInput.size() != iChroLength )
		return;
	
	// shift vdReferenceChro for the iScanShiftReference number of scans
	shiftChro( iScanShiftReference, vdReferenceChroInput );
	
	// compute noise level for treatment chro
	double dNoiseTreatmentChro = (*( min_element(vdTreatmentChroInput.begin(), 
					vdTreatmentChroInput.end() ) ));
	dNoiseTreatmentChro = dNoiseTreatmentChro - 1;

	if( dNoiseTreatmentChro < 0 )
		dNoiseTreatmentChro = 0;

	// compute noise level for reference chro
	double dNoiseReferenceChro = (*( min_element(vdReferenceChroInput.begin(), 
					vdReferenceChroInput.end() ) ));
	dNoiseReferenceChro = dNoiseReferenceChro - 1;

	if( dNoiseReferenceChro < 0 )
		dNoiseReferenceChro = 0; 

	vector<double>::iterator itr;

	// substract the noise for treatment chro
	for( itr = vdTreatmentChroInput.begin(); itr != vdTreatmentChroInput.end(); itr++ )
		*(itr) = *(itr) - dNoiseTreatmentChro;

	// substract the noise for reference chro
	for( itr = vdReferenceChroInput.begin(); itr != vdReferenceChroInput.end(); itr++ )
		*(itr) = *(itr) - dNoiseReferenceChro;

	
	// Construct covariance chromatogram
	transform( vdTreatmentChroInput.begin(), vdTreatmentChroInput.end(), vdReferenceChroInput.begin(),
			vdCovarianceChro.begin(), multiplies<double>() );

	Smoother smootherCovariance( ProRataConfig::getPeakPickerWindowSize(), ProRataConfig::getPeakPickerFOrder() );
	smootherCovariance.smoothen( vdCovarianceChro );
	
}

void PeakPicker::shiftChro( int iScanShift, vector< double > & vdChro )
{
	// no shift
	if( iScanShift == 0 )
		return;
	
	// check input validity
	if( abs( iScanShift ) > vdChro.size() - 3 || vdChro.size() < 3 )
		return;
	
	double dTempIntensity  = 0;
	
	// shift left
	if( iScanShift < 0 )
	{
		dTempIntensity = vdChro.back();
		vdChro.erase( vdChro.begin(), vdChro.begin() + abs( iScanShift ) );
		vdChro.insert( vdChro.end(), abs( iScanShift ), dTempIntensity );
	}

	// shift right
	if( iScanShift > 0 )
	{
		dTempIntensity = vdChro.front();
		vdChro.erase( ( vdChro.end() - iScanShift ) , vdChro.end() );
		vdChro.insert( vdChro.begin(), iScanShift, dTempIntensity );
	}

}

void PeakPicker::findPeak()
{

	int iOffset =  ( ProRataConfig::getPeakPickerWindowSize() - 1) / 2;
	
	// initialize iRightValley;
	for( iRightValley = ( iChroLength - 2 ) ; iRightValley > ProRataConfig::getPeakPickerWindowSize(); --iRightValley )
	{
		if ( vfRetentionTime[ iRightValley ] < fRightMS2Time )
			break;
	}
	++iRightValley;

	
	// initialize iLeftValley
	for( iLeftValley = 1 ; iLeftValley < (iChroLength - ProRataConfig::getPeakPickerWindowSize() - 1); ++iLeftValley )
	{
		if ( vfRetentionTime[ iLeftValley ] > fLeftMS2Time )
			break;
	}
	--iLeftValley;
	
	// compute iLeftValley and iRightValley
	iRightValley = findNextRightValley( iRightValley );
	iLeftValley = findNextLeftValley( iLeftValley );
	
	double dCovarianceCutOff = median( vdCovarianceChro );
	vector<double>::iterator itrBegin = vdCovarianceChro.begin();
	double dHalfPeakIntensity = ( *(max_element( itrBegin + iLeftValley, itrBegin + iRightValley + 1 )  ) - dCovarianceCutOff ) * 0.5 ;
	double dLeftValleyIntensity = vdCovarianceChro.at( iLeftValley ) - dCovarianceCutOff;
	double dRightValleyIntensity = vdCovarianceChro.at( iRightValley ) - dCovarianceCutOff;
	
	bool bMoveLeftValley = ( (dLeftValleyIntensity > dHalfPeakIntensity) && iLeftValley > iOffset );
	bool bMoveRightValley = ( (dRightValleyIntensity > dHalfPeakIntensity) && iRightValley < (iChroLength - iOffset - 1) );
	bool bMoveLeftOrRight = ( vdCovarianceChro.at(iLeftValley) > vdCovarianceChro.at(iRightValley) );

	while (	bMoveLeftValley || bMoveRightValley )
	{
		if ( bMoveLeftValley &&	( bMoveLeftOrRight || !bMoveRightValley ) )
			iLeftValley = findNextLeftValley( iLeftValley );

		if ( (!bMoveLeftValley || !bMoveLeftOrRight ) && bMoveRightValley )
			iRightValley = findNextRightValley( iRightValley );

		dHalfPeakIntensity = ( *(max_element(itrBegin + iLeftValley, itrBegin + iRightValley +1 )  ) - dCovarianceCutOff ) * 0.5 ;
		dLeftValleyIntensity = vdCovarianceChro.at( iLeftValley ) - dCovarianceCutOff;
		dRightValleyIntensity = vdCovarianceChro.at( iRightValley ) - dCovarianceCutOff;
		bMoveLeftValley = ( (dLeftValleyIntensity > dHalfPeakIntensity) && iLeftValley > iOffset ) ;
		bMoveRightValley = ( (dRightValleyIntensity > dHalfPeakIntensity) && iRightValley < (iChroLength - iOffset - 1) )	;
		bMoveLeftOrRight = ( vdCovarianceChro.at(iLeftValley) > vdCovarianceChro.at(iRightValley) );
	}

	// compute bPeakValidity 
	bool bIsPeakBelowNoise = ( *( max_element( itrBegin + iLeftValley, itrBegin + iRightValley + 1 ) ) )  < dCovarianceCutOff;
	float fPeakWidth = vfRetentionTime[ iRightValley ] - vfRetentionTime[ iLeftValley ];
//	float fChroDuration = vfRetentionTime.back() - vfRetentionTime.front() ;

	// minimal peak width = 0.1 min
	// maximal peak width = 8 min
	bool bIsPeakTooSmall = fPeakWidth < ( 0.1 );
	bool bIsPeakTooBroad = fPeakWidth > ( 8 );

	if ( bIsPeakBelowNoise || bIsPeakTooSmall || bIsPeakTooBroad )
		bPeakValidity = false;
	else
		bPeakValidity = true;


	// compute dPeakHeightPPC
	dPeakHeightPPC = ( *( max_element( itrBegin + iLeftValley, itrBegin + iRightValley + 1) ) ) -
			( ( vdCovarianceChro.at(iLeftValley) + vdCovarianceChro.at(iRightValley) )/2 );
}


double PeakPicker::median( vector<double> vdData )
{
	if( vdData.size() < 1 )
		return 0;
	
	sort( vdData.begin(), vdData.end() );
	return *( vdData.begin() + (vdData.size() / 2 ) );
}

int PeakPicker::findNextLeftValley( int iCurrentLeftValley )
{
	double dCurrentMinimumIntensity;
	double dCurrentValleyIntensity;
	int iOffset =  ( ProRataConfig::getPeakPickerWindowSize() - 1) / 2;

	bool bIsTrueValley = false;
	int iValleyMinusOffset;
	int iValleyPlusOffset;

	vector<double>::iterator itrBegin;
	itrBegin = vdCovarianceChro.begin();

	while ( (iCurrentLeftValley > iOffset ) && ( !bIsTrueValley ) )
	{
		iCurrentLeftValley--;

		dCurrentValleyIntensity = *( itrBegin + iCurrentLeftValley );

		iValleyMinusOffset = max( 0, iCurrentLeftValley - iOffset ); 
		iValleyPlusOffset = min( iChroLength -1,iCurrentLeftValley + iOffset ); 

		// Add one, because
		// min_element finds the smallest element in the range [first, last). 
		// Note the last element is not included.
		dCurrentMinimumIntensity = *( min_element(itrBegin + iValleyMinusOffset, 
					itrBegin + iValleyPlusOffset + 1 ) );

		if ( dCurrentMinimumIntensity >= dCurrentValleyIntensity )
			bIsTrueValley = true;
	}

	return iCurrentLeftValley;


}

int PeakPicker::findNextRightValley( int iCurrentRightValley )
{

	
	double dCurrentMinimumIntensity;
	double dCurrentValleyIntensity;
	int iOffset =  ( ProRataConfig::getPeakPickerWindowSize() - 1) / 2;

	bool bIsTrueValley = false;
	int iValleyMinusOffset;
	int iValleyPlusOffset;

	vector<double>::iterator itrBegin;
	itrBegin = vdCovarianceChro.begin();
	
	while ( (iCurrentRightValley < ( iChroLength - iOffset - 1 )) && 
	      		(!bIsTrueValley)  )
	{
		iCurrentRightValley++;

		dCurrentValleyIntensity = *( itrBegin + iCurrentRightValley );

		iValleyMinusOffset = max( 0, iCurrentRightValley - iOffset ); 
		iValleyPlusOffset = min( iChroLength - 1, iCurrentRightValley + iOffset ) ; 

		dCurrentMinimumIntensity = *( min_element(itrBegin + iValleyMinusOffset, 
					itrBegin + iValleyPlusOffset + 1 ) );

		if ( dCurrentMinimumIntensity >= dCurrentValleyIntensity )
			bIsTrueValley = true;
	}

	return iCurrentRightValley;
}
