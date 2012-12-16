
#include "peptideRatio.h"

PeptideRatio::PeptideRatio()
{
	dPCARatio = 1.0;
	dPCASN = 1.0;
	bValidity = true;
	iLeftValley = 0;
	iRightValley = 0;
	dReferencePeakSNR = 1.0; 
	dTreatmentPeakSNR = 1.0; 
}

PeptideRatio::~PeptideRatio()
{
	// destructor
}

bool PeptideRatio::process( const Chromatogram & inputChro )
{
	peptideChro = inputChro;

	if( !computeRatio() )
	{
		bValidity = false;
		return false;
	}

	vector< string > vsLocus;
	vector< string > vsDescription;
	if( !getLocusDescription( vsLocus, vsDescription ) )
	{
		bValidity = false;
		return false;
	}
	else
	{
		if( ( vsLocus.size() > 1 ) && ProRataConfig::getRemoveAmbiguousPeptides() )
			bValidity = false;
	}

	return true;
}

bool PeptideRatio::computeRatio()
{

	// get the chromatogram intensity for reference chro, which is the denominator by default
	if( !peptideChro.getIntensityVector( ProRataConfig::getDenominatorIsotopologue(), vdReferenceChro ) )
	{
		cout << "ERROR: cannot read Denominator Isotopologue Chromatogram! " << endl;
		return false;
	}

	// get the chromatogram intensity for treatment chro, which is the numerator by default
	if( !peptideChro.getIntensityVector( ProRataConfig::getNumeratorIsotopologue(), vdTreatmentChro ) )
	{
		cout << "ERROR: cannot read Numerator Isotopologue Chromatogram! " << endl;
		return false;
	}

	vfRetentionTime = peptideChro.getTimeVector();

	viScanVector = peptideChro.getScanVector();

	// get the vfMS2Time
	peptideChro.getMS2Time( vfMS2Time );

	// check if the chromatogram is valid
	if( !checkChromatogramValidity()  )
	{
		bValidity = false;
		return false;
	}
#ifdef DEBUG
	cout << "start peak picking " << endl;
#endif
	if( !peptidePeakPicker.process( vfMS2Time, vfRetentionTime, vdTreatmentChro, vdReferenceChro ) )
	{
		cout << "ERROR: cannot perform peak picking! " << endl;
		return false;
	}
	
	if ( peptidePeakPicker.getPeakValidity() )
	{
	
		iLeftValley = peptidePeakPicker.getLeftValley();
		iRightValley = peptidePeakPicker.getRightValley();

		if( iLeftValley >= iRightValley || iLeftValley < 0 || iRightValley >= vdReferenceChro.size())
		{
			cout << "ERROR: wrong peak boundaries!!" << endl;
			return false;
		}

		/*
		 * Calculate Peak Height and Peak Ratio
		 * Smoothing is done through Savitzky-Golay Filter
		 */
#ifdef DEBUG
	cout << "start quantification with peak height and peak area" << endl;
#endif
		if ( !quantifyWithPeak() )
		{
			cout << "Error! Could not calculate peak height and peak area ratio." << endl;
			bValidity = false;
		}

		/*
		 * Principal Component Analysis
		 */
#ifdef DEBUG
	cout << "start quantification with PCA" << endl;
#endif
		if ( !quantifyWithPCA() )
		{
		//	cout << "ERROR! Problem performing principal component analysis for Chromatogram:" << getIdentifier() << endl;
			bValidity = false;
		}

		if( dReferencePeakSNR < 1.5 && dTreatmentPeakSNR < 1.5 )
		{
			bValidity = false;
		}

		// check if this peptide's eigenvalue ratio exceed the cutoff
		// if not, set its validity to false and have it filterred out by ProteinInfo
		if( logBase2(dPCASN) < ProRataConfig::getMinLog2SNR() )
		{
			bValidity = false;
		}
			
	}
	else
	{
		bValidity = false;
	}
	
	return true;

}

/*
 * Check the validity of the data.
 * Check also the max values are non-zero.
 */

bool PeptideRatio::checkChromatogramValidity()
{	

	if( vdReferenceChro.size() < 3 || vdTreatmentChro.size() < 3 || vdTreatmentChro.size() != vdReferenceChro.size() )
	{
		bValidity = false;
		return false;
	}
	
	// Get Reference Column Max
	double dMaxReferenceValue = ( *( max_element(vdReferenceChro.begin(), vdReferenceChro.end() ) ) );
	// Get Sample Column Max
	double dMaxSampleValue = ( *( max_element(vdTreatmentChro.begin(), vdTreatmentChro.end() ) ) );

	// Check the values for zero
	if ( dMaxReferenceValue == 0 || dMaxSampleValue == 0 )
	{
	//	cout << "WARNING! chromatogram " << getIdentifier() << " has zero ion intensity" << endl;
		bValidity = false;
		return false;
	}
	return true;
}

/*
 * Perform Principal component analysis.
 * *** Explain why it is used here.
 */

bool PeptideRatio::quantifyWithPCA()
{
	// Find the width of the peak.
	int iPeakWidth = iRightValley - iLeftValley + 1;

	if( iPeakWidth < 4 )
	{
	//	cout << "ERROR: the peak width is too small! " << endl;
		bValidity = false;
		return true;
	}

	// There will be two chromatogram: treatment and reference chromatogram.
	int iNumberOfChromatograms = 2;

	int i;
	double dReferenceMean = 0;
	double dTreatmentMean = 0;

	for( i = iLeftValley; i < iRightValley + 1; ++i )
	{
		vdReferencePeak.push_back( vdReferenceChro[i] );
		dReferenceMean = dReferenceMean + vdReferenceChro[i];
		vdTreatmentPeak.push_back( vdTreatmentChro[i] );
		dTreatmentMean = dTreatmentMean + vdTreatmentChro[i]; 
//		viScanNumberPeak.push_back( peptideChro.getScanVector().at(i) );
	}

	dReferenceMean = dReferenceMean / (double) vdReferencePeak.size();
	dTreatmentMean = dTreatmentMean / (double) vdTreatmentPeak.size();

	int iNonBaselineScan = 0;
	for( i = 0; i < vdReferencePeak.size(); ++i )
	{
		if( vdReferencePeak[i] > dReferenceMean || vdTreatmentPeak[i] > dTreatmentMean )
			++iNonBaselineScan;
	}

	if( iNonBaselineScan < 2 )
	{
		bValidity = false;
		return true;
	}


	// Create an instance of PCA providing iPeakWidth as row number and 2 as column number.
	PCA chroPCA( iPeakWidth, iNumberOfChromatograms );

	/*
	for( i = 0; i < vdReferencePeak.size(); ++i )
	{
		cout << " ( " << vdReferencePeak[i] << ", " << vdTreatmentPeak[i] << " )  " << endl;
	}
	*/

	// Calculate pca on newly created vectors.
	if ( chroPCA.calculate( vdReferencePeak, vdTreatmentPeak ) != 0 )
	{
		cout << "WARNING! Calculation of Eigen values/vectors failed for Chromatogram:" << getIdentifier() << endl;
		bValidity = false;
		return false;
	}

	double dRatioUpperBound = pow( (double)2, ProRataConfig::getPCAMaxLog2Ratio() );
	double dRatioLowerBound = pow( (double)2, ProRataConfig::getPCAMinLog2Ratio() );

//	int iEigenValueCount = chroPCA.getNumberOfEigenValues();
	vector<double> vdEigenValues;
	vector<double> vdEigenVectors;
	vdEigenValues = chroPCA.getEigenValues();
	vdEigenVectors = chroPCA.getEigenVectors();

	if( vdEigenValues.size() != 2 )
	{
		cout << "ERROR: unexpected eigenvalue number for Chromatogram:" << getIdentifier() << endl;
		return false;
	}
	if( vdEigenVectors.size() != 4 )
	{
		cout << "ERROR: unexpected EigenVector number for Chromatogram:" << getIdentifier() << endl;
		return false;
	}

	if(  vdEigenValues[0] >= vdEigenValues[1] )
	{
		dFirstEigenValue = vdEigenValues[0];
		dSecondEigenValue = vdEigenValues[1];
		dxCord1 = vdEigenVectors[0];
		dyCord1 = vdEigenVectors[1];
		dxCord2 = vdEigenVectors[2];
		dyCord2 = vdEigenVectors[3];
	}
	else
	{
		dFirstEigenValue = vdEigenValues[1];
		dSecondEigenValue = vdEigenValues[0];
		dxCord1 = vdEigenVectors[2];
		dyCord1 = vdEigenVectors[3];
		dxCord2 = vdEigenVectors[0];
		dyCord2 = vdEigenVectors[1];
	}


	if ( dFirstEigenValue == 0 && dSecondEigenValue == 0 )
	{
		cout << "ERROR! both engenvalues of PCA are zero for Chromatogram:" << getIdentifier() << endl;
		bValidity = false;
		return false;
	}

	if ( dFirstEigenValue < 0 || dSecondEigenValue < 0 || dSecondEigenValue > dFirstEigenValue )
	{
		cout << "ERROR! Problematic Eigenvalues for Chromatogram:" << getIdentifier() << endl;
		cout << " first eigenvalue = " << dFirstEigenValue << endl;
		cout << " second eigenvalue = " << dSecondEigenValue << endl;
		bValidity = false;
		return false;
	}

	// use the square root the engenvalues
	dFirstEigenValue = sqrt( dFirstEigenValue );
	dSecondEigenValue = sqrt( dSecondEigenValue );
	
	
//	cout << "PCA  X1 = " << dxCord1 << " Y1 = " << dyCord1 << " X2 = " << dxCord2 << " Y2 =  "  <<  dyCord2 << endl;
//	cout << "dFirstEigenValue = " << dFirstEigenValue << "  dSecondEigenValue = "  << dSecondEigenValue << endl;

	/*
	if ( (dxCord1 * dxCord2 * dyCord1 * dyCord2) == 0 )
	{
		if ( dFirstEigenValue >= dSecondEigenValue )
		{
			if ( (dxCord1 * dyCord1) > 0 )
			{
				dPCARatio = dRatioUpperBound; 
			}

			if ( (dxCord2 * dyCord2) > 0 )
			{
				dPCARatio = dRatioLowerBound; 
			}

			if( dSecondEigenValue != 0 )
			{
				dPCASN = dFirstEigenValue / dSecondEigenValue;
			}
			else
			{
				dPCASN = 1024;
			}
		}
		else
		{
			if ( (dxCord2 * dyCord2) > 0 )
			{
				dPCARatio = dRatioUpperBound; 
			}

			if ( (dxCord1 * dyCord1) > 0 )
			{
				dPCARatio = dRatioLowerBound;
			}
			
			if( dFirstEigenValue != 0 )
			{
				dPCASN = dSecondEigenValue / dFirstEigenValue;;
			}
			else
			{
				dPCASN = 1024;
			}
						
		}
		
	}
	else
	{
		
		if ( (dxCord1 * dyCord1) > 0 )
		{
			dPCARatio = dyCord1 / dxCord1;
			if( dSecondEigenValue != 0 )
			{
				dPCASN = dFirstEigenValue / dSecondEigenValue;
			}
			else
			{
				dPCASN = 1024;
			}
		}

		if ( (dxCord2 * dyCord2) > 0 )
		{
			dPCARatio = dyCord2 / dxCord2;
			if( dFirstEigenValue != 0 )
			{
				dPCASN = dSecondEigenValue / dFirstEigenValue;
			}
			else
		vdPeptideLog2RatioInput[i]	{
				dPCASN = 1024;
			}
		}

		if ( dPCARatio > dRatioUpperBound || dPCARatio < -1024  )
		{
			dPCARatio = dRatioUpperBound;
		}

		if ( ( -1.0/1024.0) < dPCARatio && dPCARatio < dRatioLowerBound )
		{
			dPCARatio = dRatioLowerBound;
		}

	}
	*/
	if ( dxCord1 == 0 )
	{
		dPCARatio = dRatioUpperBound;
	}
	else
	{
		dPCARatio = dyCord1 / dxCord1;

		if ( dPCARatio > dRatioUpperBound || dPCARatio < -1024  )
		{
			dPCARatio = dRatioUpperBound;
		}

		if ( ( -1.0/1024.0) < dPCARatio && dPCARatio < dRatioLowerBound )
		{
			dPCARatio = dRatioLowerBound;
		}

	}

	if( dPCARatio < 0 )
	{
		dPCARatio = 1;
		dPCASN = 1;

	}
	else
	{
		if( dSecondEigenValue != 0 )
		{
			dPCASN = dFirstEigenValue / dSecondEigenValue;
		}
		else
		{
			dPCASN = 1024;
		}
	}
	
			
	return true;
}


/*
 * Calculate Peak height ratio and peak area ratio.
 * *** Explain why they are significant.
 */

bool PeptideRatio::quantifyWithPeak()
{

	Smoother peakSmoother( ProRataConfig::getPeakHeightWindowSize(), ProRataConfig::getPeakHeightFOrder() );
	
	vector< double > vdTreatmentChroCopy = vdTreatmentChro;
	vector< double > vdReferenceChroCopy = vdReferenceChro;

	// Smoothen vdTreatmentChroCopy.
	peakSmoother.smoothen( vdTreatmentChroCopy );

	// Smoothen vdTreatmentChroCopy
	peakSmoother.smoothen( vdReferenceChroCopy );

	// compute the peak SNR
	dReferencePeakSNR = computePeakSNR( vdReferenceChroCopy );
	dTreatmentPeakSNR = computePeakSNR( vdTreatmentChroCopy );


	// Find the size of the valley.
	int iPeakWidth = iRightValley - iLeftValley + 1;
	
	int i;

	// get the peak intensities within a peak
	vector<double> vdTreatmentIntensity;
	vector<double> vdReferenceIntensity;	
	for( i = iLeftValley; i < iRightValley + 1; ++i )
	{
		vdReferenceIntensity.push_back( vdReferenceChroCopy[i] );
		vdTreatmentIntensity.push_back( vdTreatmentChroCopy[i] );
	}

	// Create two vectors to store semi-synthetic base values.
	vector<double> vdTreatmentBaseline( iPeakWidth );
	vector<double> vdReferenceBaseline( iPeakWidth );
	
	double dFirstValue = vdTreatmentIntensity.front();
	double dLastValue = vdTreatmentIntensity.back();
	double dTreatmentDelta = ( dLastValue - dFirstValue ) / ( iPeakWidth - 1 );
	
	vdTreatmentBaseline.at( 0 ) = dFirstValue;
	for( i = 1; i < vdTreatmentIntensity.size(); i++ )
	{
		vdTreatmentBaseline.at( i ) = vdTreatmentBaseline.at( i - 1 ) + dTreatmentDelta;
	}


	dFirstValue = vdReferenceIntensity.front();
	dLastValue = vdReferenceIntensity.back();
	double dReferenceDelta =  ( dLastValue - dFirstValue ) / ( iPeakWidth - 1 );

	vdReferenceBaseline.at( 0 ) = dFirstValue;
	for( i = 1; i < vdReferenceIntensity.size() ; i++ )
	{
		vdReferenceBaseline.at( i ) = vdReferenceBaseline.at( i - 1 ) + dReferenceDelta;
	}


	for( i = 0; i < vdTreatmentIntensity.size(); ++i )
	{
		vdTreatmentIntensity[i] = vdTreatmentIntensity[i] - vdTreatmentBaseline[i];
		if( vdTreatmentIntensity[i] < 0 )
			vdTreatmentIntensity[i] = 0;
		vdReferenceIntensity[i] = vdReferenceIntensity[i] - vdReferenceBaseline[i];
		if( vdReferenceIntensity[i] < 0 )
			vdReferenceIntensity[i] = 0;
	}

	double dTreatmentPeakArea = accumulate( vdTreatmentIntensity.begin(), vdTreatmentIntensity.end(), 0.0 );
	double dTreatmentPeakHeight = ( *( max_element( vdTreatmentIntensity.begin(), vdTreatmentIntensity.end() ) ) );

	double dReferencePeakArea = accumulate( vdReferenceIntensity.begin(), vdReferenceIntensity.end(), 0.0 );
	double dReferencePeakHeight = ( *( max_element( vdReferenceIntensity.begin(), vdReferenceIntensity.end() ) ) );

	if ( dTreatmentPeakArea > 0 && dReferencePeakArea > 0 )
	{
		dPeakAreaRatio = dTreatmentPeakArea / dReferencePeakArea ;
	}
	else if( dTreatmentPeakArea <= 0 )
	{
		dPeakAreaRatio = ProRataConfig::getPCAMinLog2Ratio();
	}
	else 
	{
		dPeakAreaRatio = ProRataConfig::getPCAMaxLog2Ratio();
	}

	if ( dTreatmentPeakHeight > 0 && dReferencePeakHeight > 0 )
	{
		dPeakHeightRatio = dTreatmentPeakHeight / dReferencePeakHeight ;
	}
	else if( dTreatmentPeakHeight <= 0 )
	{
		dPeakHeightRatio = ProRataConfig::getPCAMinLog2Ratio();
	}
	else
	{
		dPeakHeightRatio = ProRataConfig::getPCAMaxLog2Ratio();
	}

	return true;
}



double PeptideRatio::computePeakSNR( vector< double > vdChro )
{
	double dPeakSNR = 1;
	double dPeak = 1;
	double dNoise = 1;

	
	dPeak = ( *( max_element(vdChro.begin() + iLeftValley, vdChro.begin() + iRightValley + 1) ) );
	vdChro.erase( (vdChro.begin() + iLeftValley), (vdChro.begin() + iRightValley + 1) );

	if( vdChro.size() < 3 )
	       return dPeakSNR;	
/*			
	int i;
 	double dBaseline = 0;
	for( i = 0; i < vdChro.size(); ++i )
	{
		dBaseline = dBaseline + vdChro[i];
	}
	dBaseline = dBaseline / (double)vdChro.size();

	dPeak = dPeak - dBaseline;

	double dSD = 0;
	for( i = 0; i < vdChro.size(); ++i )
	{
		dSD = dSD + ( vdChro[i] - dBaseline ) * ( vdChro[i] - dBaseline ); 
	}

	dSD = dSD / (double)( vdChro.size() - 1 );

	dSD = sqrt( dSD );

	dNoise = dSD;
*/
	
	sort( vdChro.begin(), vdChro.end() );
//	double dNoiseMinimum = *( vdChro.begin() + (int)(vdChro.size() * 0.1 ) ); 
//	double dNoiseMaximum = *( vdChro.begin() + (int)(vdChro.size() * 0.9 ) ); 
//	double dNoise = ( dNoiseMaximum - dNoiseMinimum );
	double dNoiseMedian = *( vdChro.begin() + (int)(vdChro.size() *0.5 ) );
	double dNoiseMinimum = *( vdChro.begin() + (int)(vdChro.size() * 0.1 ) ); 
	dNoise = ( dNoiseMedian - dNoiseMinimum );
	dPeak = dPeak - dNoiseMedian;
	
	if( dPeak > 0 && dNoise > 0 )
		dPeakSNR = dPeak / dNoise;
	else if( dPeak > 0 && dNoise <= 0 )
		dPeakSNR = 1024;
	else
		dPeakSNR = 1;

	return dPeakSNR;
}

bool PeptideRatio::getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription )
{
	vector< string > vsLocusTemp;
	vector< string > vsDescriptionTemp;

	bool bFlag = peptideChro.getLocusDescription( vsLocusTemp, vsDescriptionTemp );

	vsLocus = vsLocusTemp;
	vsDescription = vsDescriptionTemp;

	return bFlag;
}


bool PeptideRatio::getReferenceMZrange( vector< float > & vfLowerMZ, vector< float > & vfUpperMZ )
{
	vector< float > vfLowerMZTemp;
	vector< float > vfUpperMZTemp;
	// get the m/z range for reference chro, which is the denominator by default
	if( !peptideChro.getMZwindows( ProRataConfig::getDenominatorIsotopologue(), vfLowerMZTemp, vfUpperMZTemp ) )
	{
		cout << "ERROR: cannot read Denominator Isotopologue Chromatogram! " << endl;
		return false;
	}
	vfLowerMZ.clear();
	vfUpperMZ.clear();

	vfLowerMZ = vfLowerMZTemp;
	vfUpperMZ = vfUpperMZTemp;

	return true;

}

bool PeptideRatio::getTreatmentMZrange( vector< float > & vfLowerMZ, vector< float > & vfUpperMZ )
{
	
	vector< float > vfLowerMZTemp;
	vector< float > vfUpperMZTemp;
	// get the m/z range for treatment chro, which is the numerator by default
	if( !peptideChro.getMZwindows( ProRataConfig::getNumeratorIsotopologue(), vfLowerMZTemp, vfUpperMZTemp ) )
	{
		cout << "ERROR: cannot read Numerator Isotopologue Chromatogram! " << endl;
		return false;
	}
	vfLowerMZ.clear();
	vfUpperMZ.clear();

	vfLowerMZ = vfLowerMZTemp;
	vfUpperMZ = vfUpperMZTemp;

	return true;
}

void PeptideRatio::getEngenvectors( double & dPC1x, double & dPC1y, double & dPC2x, double & dPC2y )
{
	dPC1x = dxCord1; 
	dPC1y = dyCord1;	
	dPC2x = dxCord2;
	dPC2y = dyCord2;
	return;
}

void PeptideRatio::getEngenvalues( double & dPC1EV, double & dPC2EV )
{
	dPC1EV = dFirstEigenValue;
	dPC2EV = dSecondEigenValue;
	return;
}

vector< unsigned long int > PeptideRatio::getScanNumberPeak()
{
	vector< unsigned long int > viScanNumberPeak;
	int i;
	for( i = iLeftValley; i < iRightValley + 1; ++i )
	{
		viScanNumberPeak.push_back( viScanVector[i] );
	}

	return viScanNumberPeak;
}

unsigned long int PeptideRatio::getNextFullScan( unsigned long int iCurrentFullScan, bool bMoveForward )
{
	if( iCurrentFullScan == 0 )
		return viScanVector[iLeftValley];

	if( iCurrentFullScan < viScanVector.front() )
		return viScanVector.front();

	if( iCurrentFullScan == viScanVector.front() && !bMoveForward )
		return viScanVector.back();

	if( iCurrentFullScan > viScanVector.back() )
		return viScanVector.back();

	if( iCurrentFullScan == viScanVector.back() && bMoveForward )
		return viScanVector.front();

	unsigned int iScanSize = viScanVector.size();
	if( bMoveForward )
	{
		for( unsigned int i = 0; i < iScanSize; ++i )
		{
			if( viScanVector[i] <= iCurrentFullScan && iCurrentFullScan < viScanVector[i+1] )
			{
				return viScanVector[i+1];
			}
		}
	}
	else
	{
		for( unsigned int i = 0; i < iScanSize; ++i )
		{
			if( viScanVector[i] < iCurrentFullScan && iCurrentFullScan <= viScanVector[i+1] )
			{
				return viScanVector[i];
			}
		}

	}

	return viScanVector.front();

}

float PeptideRatio::getFullScanTime( unsigned long int iScan )
{
	unsigned int iScanSize = viScanVector.size();
	for( unsigned int i = 0; i < iScanSize; ++i )
	{
		if( viScanVector[i] >= iScan )
		{
			return vfRetentionTime[i];
		}
	}
	return vfRetentionTime.back();
}
