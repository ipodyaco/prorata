#ifndef PROTEINREPLICATE_H
#define PROTEINREPLICATE_H

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <iomanip>
#include <fstream>

#include "proteinInfo.h"
#include "proteomeInfo.h"
#include "proteinRatio.h"
#include "directoryStructure.h"
#include "proRataConfig.h"
#include "tinyxml.h"

using namespace std;

class ProteinReplicate
{
	
	public:
		ProteinReplicate();
		~ProteinReplicate();

		// calculate results
		bool runQuantification( string sInputLocus );
	       	void eraseProfileLikelihoodCurve();
		
		// get results
		bool getProfileLikelihoodCurve(vector<double> & vdLog2RatioRef, vector<double> & vdLnLikelihoodRef );
		void getEstimates(bool & bValidityRef, double & dLog2RatioRef, double & dLowerLimitCIRef, double & dUpperLimitCIRef );
		void getEstimates(double & dLowerLimitCIRef, double & dUpperLimitCIRef );
		int getPeptideNumberWithinCI( double dLowerLimitCI, double dUpperLimitCI );
		int getQuantifiedPeptides();
		bool isValid();	
		string getLocus();

		// these two functions are called only in projectInfo to setupTemplateProteinCombined
		void setName( string sNameConfig );
		string getName();
		void setProteomeInfo( ProteomeInfo * pProteomeInfoConfig );

		// get CONFIG
		void getLog2RatioSetting(double & dMLEMinLog2RatioRef, double & dMLEMaxLog2RatioRef, double & dLog2RatioDiscretizationRef);

	private:

		// read from ProRataCombine.XML
		string sName;
		ProteomeInfo * pProteomeInfo;

		// variables from the PROTEIN_QUANTIFICATION element
		// copied from proRataConfig.h
		bool bRemoveAmbiguousPeptide;
		int iMinPeptideNumber;
		double dMaxCIwidth;
		double dMinLog2SNR;
		double dMaxLog2SNR;
		double dMLEMinLog2Ratio;
		double dMLEMaxLog2Ratio;
		double dLog2RatioDiscretization;
		double dSDSlope;
		double dSDIntercept;
		double dMeanSlope;
		double dMeanIntercept;
		double dSmoothingProbSpace;
		double dLnLikelihoodCutoffOffset;
		
		// locus and corresponding pProteinInfo of the protein
		string sLocus;
		ProteinInfo * pProteinInfo;

		// estimaiton results
		bool bValidity;
		double dLog2Ratio;
		double dLowerLimitCI;
		double dUpperLimitCI;
		int iQuantifiedPeptides;

		// profile likelihood curve
		vector<double> vdLog2Ratio;
		vector<double> vdLnLikelihood;

};

#endif //PROTEINREPLICATE_H
