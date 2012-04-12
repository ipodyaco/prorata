
#include "peptideLabelFree.h"

PeptideLabelFree::PeptideLabelFree()
{
	bValidity = true;
	iLeftValley = 0;
	iRightValley = 0;
	dReferencePeakSNR = 1.0; 
	dPeakHeight = 0;
	dPeakArea = 0;	
}

PeptideLabelFree::~PeptideLabelFree()
{
	// destructor
}

bool PeptideLabelFree::process( const Chromatogram & inputChro )
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

bool PeptideLabelFree::computeRatio()
{
	// get the chromatogram
	vector<SIC> vSIC = peptideChro.getAllSIC();
	if( vSIC.size() != 1 )
	{
		cout << "ERROR: There should be one and only one isotopologue for label-free analysis." << endl;
		return false;
	}
	vdReferenceChro = vSIC[0].vdIntensity;
	sChroName = vSIC[0].sName;

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

	if( !peptidePeakPicker.process( vfMS2Time, vfRetentionTime, vdReferenceChro, vdReferenceChro ) )
	{
		cout << "ERROR: cannot perform peak picking! " << endl;
		return false;
	}
	
	if ( !peptidePeakPicker.getPeakValidity() )
	{
		bValidity = false;		
		return true;
	}
	
	iLeftValley = peptidePeakPicker.getLeftValley();
	iRightValley = peptidePeakPicker.getRightValley();

	if( iLeftValley >= iRightValley || iLeftValley < 0 || iRightValley >= (int)vdReferenceChro.size())
	{
		cout << "ERROR: wrong peak boundaries!!" << endl;
		return false;
	}

	/*
	 * Calculate Peak Height and Peak Ratio
	 * Smoothing is done through Savitzky-Golay Filter
	 */
	if ( !quantifyWithPeak() )
	{
		cout << "Error! Could not calculate peak height and peak area ratio." << endl;
		bValidity = false;
	}
	
	return true;

}

/*
 * Check the validity of the data.
 * Check also the max values are non-zero.
 */

bool PeptideLabelFree::checkChromatogramValidity()
{	
	if( vdReferenceChro.size() < 3 )
	{
		bValidity = false;
		return false;
	}
	
	// Get Reference Column Max
	double dMaxReferenceValue = ( *( max_element(vdReferenceChro.begin(), vdReferenceChro.end() ) ) );

	// Check the values for zero
	if ( dMaxReferenceValue == 0 )
	{
	//	cout << "WARNING! chromatogram " << getIdentifier() << " has zero ion intensity" << endl;
		bValidity = false;
		return false;
	}
	return true;
}

/*
 * Calculate Peak height ratio and peak area ratio.
 * *** Explain why they are significant.
 */

bool PeptideLabelFree::quantifyWithPeak()
{

	Smoother peakSmoother( ProRataConfig::getPeakHeightWindowSize(), ProRataConfig::getPeakHeightFOrder() );
	
	vector< double > vdReferenceChroCopy = vdReferenceChro;

	// Smoothen vdReferenceChroCopy
	peakSmoother.smoothen( vdReferenceChroCopy );

	// compute the peak SNR
	dReferencePeakSNR = computePeakSNR( vdReferenceChroCopy );

	// Find the size of the valley.
	int iPeakWidth = iRightValley - iLeftValley + 1;
	
	unsigned int i;

	// get the peak intensities within a peak
	vector<double> vdReferenceIntensity;	
	for( i = iLeftValley; i < iRightValley + 1; ++i )
	{
		vdReferenceIntensity.push_back( vdReferenceChroCopy[i] );
	}

	// Create two vectors to store semi-synthetic base values.
	vector<double> vdReferenceBaseline( iPeakWidth );
	
	double dFirstValue = vdReferenceIntensity.front();
	double dLastValue = vdReferenceIntensity.back();
	double dReferenceDelta =  ( dLastValue - dFirstValue ) / ( iPeakWidth - 1 );

	vdReferenceBaseline.at( 0 ) = dFirstValue;
	for( i = 1; i < vdReferenceIntensity.size() ; i++ )
	{
		vdReferenceBaseline.at( i ) = vdReferenceBaseline.at( i - 1 ) + dReferenceDelta;
	}


	for( i = 0; i < vdReferenceIntensity.size(); ++i )
	{
		vdReferenceIntensity[i] = vdReferenceIntensity[i] - vdReferenceBaseline[i];
		if( vdReferenceIntensity[i] < 0 )
			vdReferenceIntensity[i] = 0;
	}

	dPeakArea = accumulate( vdReferenceIntensity.begin(), vdReferenceIntensity.end(), 0.0 );
	dPeakHeight = ( *( max_element( vdReferenceIntensity.begin(), vdReferenceIntensity.end() ) ) );

	return true;
}



double PeptideLabelFree::computePeakSNR( vector< double > vdChro )
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

bool PeptideLabelFree::getLocusDescription( vector< string > & vsLocus, vector< string > & vsDescription )
{
	vector< string > vsLocusTemp;
	vector< string > vsDescriptionTemp;

	bool bFlag = peptideChro.getLocusDescription( vsLocusTemp, vsDescriptionTemp );

	vsLocus = vsLocusTemp;
	vsDescription = vsDescriptionTemp;

	return bFlag;
}

vector< unsigned long int > PeptideLabelFree::getScanNumberPeak()
{
	vector< unsigned long int > viScanNumberPeak;
	int i;
	for( i = iLeftValley; i < iRightValley + 1; ++i )
	{
		viScanNumberPeak.push_back( viScanVector[i] );
	}

	return viScanNumberPeak;
}

