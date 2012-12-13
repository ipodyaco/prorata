
#ifndef PRORATACONFIG_H
#define PRORATACONFIG_H

#include "tinyxml.h"
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;
typedef map< string, string, less<string> > residueMap;

class ProRataConfig
{
	public:
		/*
		 * Sets up sessionwide configuration for ProRata
		 * the configurations are loaded in to memory as static variables
		 */

		static bool setFilename( const string & sConfigFileName );

		static bool setWorkingDirectory( const string & sDirectoryName );

		// whether or not to write chro files
		static void setWriteChro(bool input){bIfWriteChro = input;};
		static bool getWriteChro(){return bIfWriteChro;};
		
		static string getWorkingDirectory()
		{ return sWorkingDirectory; }
		
		static char getSeparator();

		// If there is only one isotopologue, this is a label-free run and return true
		static bool getIsLabelFree() {return bIsLabelFree;};

		/*
		 * get the version number of ProRata
		 */
		
		static string getProRataVersion() { return "3.0"; }
		
		/*
		 * the get functions for SIC extraction 
		 * <SIC_EXTRACTION>
		 */ 
		
		// retrieve <RETENTION_TIME_INTERVAL> <MINUTES_BEFORE_MS2>
		static float	getMinutesBeforeMS2()			{ return fMinutesBeforeMS2; }
		static void	setMinutesBeforeMS2(float input)	{ fMinutesBeforeMS2 = input; }
		
		// retrieve <RETENTION_TIME_INTERVAL> <MINUTES_AFTER_MS2>
		static float	getMinutesAfterMS2()			{ return fMinutesAfterMS2; }
		static void	setMinutesAfterMS2(float input)		{ fMinutesAfterMS2 = input; }
		
		// retrieve <RETENTION_TIME_INTERVAL> <MINUTES_BETWEEN_DUPLICATE_MS2>
		static float	getMinutesBetweenMS2()			{ return fMinutesBetweenMS2; }
		static void	setMinutesBetweenMS2(float input)	{ fMinutesBetweenMS2 = input; }
		
		// retrieve <MASS_TO_CHARGE_INTERVAL> <PLUS_MZ_ERROR>
		static double	getPlusMZerror()			{ return dPlusMZerror; }
		static void	setPlusMZerror(double input)		{ dPlusMZerror = input; }
		
		// retrieve <MASS_TO_CHARGE_INTERVAL> <MINUS_MZ_ERROR>
		static double	getMinusMZerror()			{ return dMinusMZerror; }
		static void	setMinusMZerror(double input)		{ dMinusMZerror = input; }
		
		// retrieve <MASS_TO_CHARGE_INTERVAL> <ISOTOPIC_ENVELOP_CUTOFF>
		static double	getIsotopicEnvelopCutoff()		{ return dIsotopicEnvelopCutoff; }
		static void	setIsotopicEnvelopCutoff(double input)	{ dIsotopicEnvelopCutoff = input; }
		
		// retrieve <ATOM_ISOTOPIC_COMPOSITION>
		// the input character is the atom name CHONPS
		static bool getAtomIsotopicComposition(
				char cAtom, 
				vector<double> & vdAtomicMass,  
				vector<double> & vdNaturalComposition,
				vector<double> & vdEnrichedComposition);
		
		// retrieve <RESIDUE_ATOMIC_COMPOSITION>
		static bool getResidueAtomicComposition( residueMap & mIsotopologue );
		
		/*
		 * the get functions for peptide quantification
		 * <PEPTIDE_QUANTIFICATION>
		 */
		
		// retrieve <PEAK_DETECTION> <CHROMATOGRAM_SMOOTHING> <WINDOW_SIZE>
		static int	getPeakPickerWindowSize()		{ return iPeakDetectionWindowSize; }
		static void	setPeakPickerWindowSize(int input)	{ iPeakDetectionWindowSize = input; }

		// retrieve <PEAK_DETECTION> <CHROMATOGRAM_SMOOTHING> <ORDER>
		static int	getPeakPickerFOrder()			{ return iPeakDetectionFOrder; }
		static void	setPeakPickerFOrder(int input)		{ iPeakDetectionFOrder = input; }

		// the parameters for Sav-Gol smoothing in peak height/are calculation
		// they will be removed
		static int	getPeakHeightWindowSize()		{ return iPeakDetectionWindowSize; }
		static void	setPeakHeightWindowSize(int input)	{ iPeakDetectionWindowSize = input; }
		
		static int	getPeakHeightFOrder()			{ return iPeakDetectionFOrder; }
		static void	setPeakHeightFOrder(int input)		{ iPeakDetectionFOrder = input; }
		
		// retrieve <PEAK_DETECTION> <PEAK_SHIFT> <LEFT>
		static float	getLeftPeakShift()			{ return fLeftPeakShift; }
		static void	setLeftPeakShift(float input)		{ fLeftPeakShift = input; }

		// retrieve <PEAK_DETECTION> <PEAK_SHIFT> <RIGHT>
		static float	getRightPeakShift()			{ return fRightPeakShift; }
		static void	setRightPeakShift(float input)		{ fRightPeakShift = input; }
		
		// retrieve <ABUNDANCE_RATIO> <NUMERATOR_ISOTOPOLOGUE>
		static string	getNumeratorIsotopologue()		{ return sNumeratorIsotopologue; }
		static void	setNumeratorIsotopologue(string input)	{ sNumeratorIsotopologue = input; }
		
		// retrieve <ABUNDANCE_RATIO> <DENOMINATOR_ISOTOPOLOGUE>
		static string	getDenominatorIsotopologue()		{ return sDenominatorIsotopologue; }
		static void	setDenominatorIsotopologue(string input) { sDenominatorIsotopologue = input; }
		
		// retrieve <LOG2_RATIO> <MINIMUM>
		static double	getPCAMinLog2Ratio()			{ return dPCAMinLog2Ratio; }
		static void	setPCAMinLog2Ratio(double input)	{ dPCAMinLog2Ratio = input; }
		
		// retrieve <LOG2_RATIO> <MAXIMUM>
		static double	getPCAMaxLog2Ratio()			{ return dPCAMaxLog2Ratio; }
		static void	setPCAMaxLog2Ratio(double input)	{ dPCAMaxLog2Ratio = input; }
	
		/* 
		 * the get functions for protein quantification
		 * <PROTEIN_QUANTIFICATION> 
		 */
		
		// retrieve <REMOVE_AMBIGUOUS_PEPTIDES>	
		static bool	getRemoveAmbiguousPeptides()		{ return bRemoveAmbiguousPeptide; }
		static void	setRemoveAmbiguousPeptides(bool input)	{ bRemoveAmbiguousPeptide = input; }

		// retrieve <MIN_PEPTIDE_NUMBER>
		static int	getMinPeptideNumber()			{ return iMinPeptideNumber; }
		static void	setMinPeptideNumber(int input)		{ iMinPeptideNumber = input; }

		// retrieve <MAX_CI_WIDTH>
		static double	getMaxCIwidth()				{ return dMaxCIwidth; }
		static void	setMaxCIwidth(double input)		{ dMaxCIwidth = input; }
		
		// retrieve <LOG2_SNR> <MINIMUM>
		static double	getMinLog2SNR()				{ return dMinLog2SNR; }
		static void	setMinLog2SNR(double input)		{ dMinLog2SNR = input; }

		// retrieve <LOG2_SNR> <MAXIMUM>
		static double	getMaxLog2SNR()				{ return dMaxLog2SNR; }
		static void	setMaxLog2SNR(double input)		{ dMaxLog2SNR = input; }

		// retrieve <LOG2_RATIO> <MINIMUM>
		static double	getMLEMinLog2Ratio()			{ return dMLEMinLog2Ratio; }
		static void	setMLEMinLog2Ratio(double input)	{ dMLEMinLog2Ratio = input; }
		
		// retrieve <LOG2_RATIO> <MAXIMUM>
		static double	getMLEMaxLog2Ratio()			{ return dMLEMaxLog2Ratio; }
		static void	setMLEMaxLog2Ratio(double input)	{ dMLEMaxLog2Ratio = input; }

		// retrieve <LOG2_RATIO_DISCRETIZATION>
		static double	getLog2RatioDiscretization()		{ return dLog2RatioDiscretization; }
		static void	setLog2RatioDiscretization(double input)	{ dLog2RatioDiscretization = input; }

		// retrieve <STANDARD_DEVIATION> <SLOPE>
		static double	getSDSlope()				{ return dSDSlope; }
		static void	setSDSlope(double input)		{ dSDSlope = input; }

		// retrieve <STANDARD_DEVIATION> <INTERCEPT>
		static double	getSDIntercept()			{ return dSDIntercept; }
		static void	setSDIntercept(double input)		{ dSDIntercept = input; }
		
		// retrieve <MEAN> <SLOPE>
		static double	getMeanSlope()				{ return dMeanSlope; }
		static void	setMeanSlope(double input)		{ dMeanSlope = input; }
		
		// retrieve <MEAN> <INTERCEPT>
		static double	getMeanIntercept()			{ return dMeanIntercept; }
		static void	setMeanIntercept(double input)		{ dMeanIntercept = input; }
			
		// retrieve <SMOOTHING_PROBABILITY_SPACE>
		static double	getSmoothingProbilitySpace()		{ return dSmoothingProbSpace; }
		static void	setSmoothingProbilitySpace(double input)	{ dSmoothingProbSpace = input; }

		// retrieve <LN_LIKELIHOOD_CUTOFF_OFFSET>
		static double	getLnLikelihoodCutoffOffset()		{ return dLnLikelihoodCutoffOffset; }
		static void	setLnLikelihoodCutoffOffset(double input)	{ dLnLikelihoodCutoffOffset = input; }

		// write the CONFIG to a file with the tab depth of iTabDepth
		// the booleans for each element controls which element to be written to the file
		static bool writeConfigXML( string sOutputFile, int iTabDepth, 
				bool bSicExtract, bool bPeptideQuan, bool bProteinQuan );

	protected:
		ProRataConfig();

	private:

		static ProRataConfig* proRataConfigSingleton;

		void setIsLabelFree();

		// the filename of the configuration file
		static string sFilename;

		// the working directory
		static string sWorkingDirectory;

		// whether or not to write chro files
		static bool bIfWriteChro;

		// whether this is a label-free analysis or not
		static bool bIsLabelFree;

		// load all static parameters into the memory
		// when openning the file 
		void getParameters( TiXmlDocument & );

		// retrieve the value as string from an element
		// the element is located by giving its hierarchial path
		static string getValue( TiXmlDocument &, const vector<string>& );

		static void replaceDelimitor( string & sLine, char cOldDelimitor, char cNewDelimitor );

		// variables from the SIC_EXTRACTION element
		static float fMinutesBeforeMS2;
		static float fMinutesAfterMS2;
		static float fMinutesBetweenMS2;
		static double dPlusMZerror;
		static double dMinusMZerror;
		static double dIsotopicEnvelopCutoff;

		// variables from the PEPTIDE_QUANTIFICATION element
		static int iPeakDetectionWindowSize;
		static int iPeakDetectionFOrder;
		
		static float fLeftPeakShift;
		static float fRightPeakShift;

		static string sNumeratorIsotopologue;
		static string sDenominatorIsotopologue;

		static double dPCAMinLog2Ratio;
		static double dPCAMaxLog2Ratio;

		// variables from the PROTEIN_QUANTIFICATION element
		
		static bool bRemoveAmbiguousPeptide;

		static int iMinPeptideNumber;

		static double dMaxCIwidth;
		
		static double dMinLog2SNR;
		static double dMaxLog2SNR;

		static double dMLEMinLog2Ratio;
		static double dMLEMaxLog2Ratio;

		static double dLog2RatioDiscretization;

		static double dSDSlope;
		static double dSDIntercept;
		
		static double dMeanSlope;
		static double dMeanIntercept;

		static double dSmoothingProbSpace;

		static double dLnLikelihoodCutoffOffset;


};
#endif //PRORATACONFIG_H
