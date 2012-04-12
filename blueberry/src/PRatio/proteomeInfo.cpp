
#include "proteomeInfo.h"

ProteomeInfo::ProteomeInfo()
{
	// constructor
}

ProteomeInfo::~ProteomeInfo()
{
	for( int i = 0; i < vpProteinInfo.size(); ++i )
		delete vpProteinInfo[i];

	for( int j = 0; j < vpPeptideInfo.size(); ++j )
		delete vpPeptideInfo[j];
}

bool ProteomeInfo::processPeptidesXIC(string sIDFilename)
{
	SICinfo sicInfo;
	sicInfo.setFilename( sIDFilename );
	sicInfo.process( vpPeptideInfo);

	for( int i = 0; i < vpPeptideInfo.size(); i++ )
	{
		addPeptideInfo( vpPeptideInfo[i] );
	}
	return true;
}

void ProteomeInfo::addPeptideInfo( PeptideInfo * pCurrentPeptideInfo )
{
	vector< string > vsLocus = pCurrentPeptideInfo->getLocus();
	vector< string > vsDescription = pCurrentPeptideInfo->getDescription();
	bool bIsNewLocus = true;
	int iLocusNumber = vsLocus.size();
	int i;
	int j;
	for( i = 0; i < iLocusNumber; ++i )
	{
		bIsNewLocus = true;
		for( j = 0; j < vpProteinInfo.size(); ++j )
		{
			if( vsLocus[i] == vpProteinInfo[j]->getLocus() )
			{
				bIsNewLocus = false;
				vpProteinInfo[j]->addPeptideInfo( pCurrentPeptideInfo );
			}
		}
		if( bIsNewLocus )
		{
			ProteinInfo * pProteinInfo = new ProteinInfo;
			pProteinInfo->setLocus( vsLocus[i] );
			pProteinInfo->setDescription( vsDescription[i] );
			pProteinInfo->addPeptideInfo( pCurrentPeptideInfo );
			vpProteinInfo.push_back( pProteinInfo );
		}
	}

}

bool ProteomeInfo::processProteins()
{
	cout << " Calculating Protein Abundance Ratios ... " << endl;
	unsigned int i;
	int iCount = 0;
	if(!ProRataConfig::getIsLabelFree())
	{
		// calculate abundance ratios for metabolic labeling
		for( i = 0; i < vpProteinInfo.size(); ++i )
		{
		//	cout << "processing " << vpProteinInfo[i]->getLocus() << endl;
			ProteinRatio currentProteinRatio;
			vpProteinInfo[i]->setProteinRatio( & currentProteinRatio );
			// count how many invalid ProteinInfo
			if( !vpProteinInfo[i]->getValidity() )
			{
				++iCount;
			}
		}
	}
	else
	{
		// perform label-free quantification
		for( i = 0; i < vpProteinInfo.size(); ++i )
		{
		//	cout << "processing " << vpProteinInfo[i]->getLocus() << endl;
			vpProteinInfo[i]->computeLabelFree();
			// count how many invalid ProteinInfo
			if( !vpProteinInfo[i]->getValidity() )
			{
				++iCount;
			}
		}

	}

	// erase those invalid ProteinInfo
	// change later to use erase_if
	vector< ProteinInfo * >::iterator iter;
	for( i = 1; i < (iCount+1); ++i )
	{
		for( iter = vpProteinInfo.begin(); iter != vpProteinInfo.end(); iter++ )
		{
			if( !(*iter)->getValidity() )
			{
				delete (*iter);
				vpProteinInfo.erase( iter );
				break;
			}
		}
	}
	
	return true;
}

bool ProteomeInfo::getProteinRatio( ProteinInfo * queryProteinInfo, ProteinRatio * queryProteinRatio )
{
	return queryProteinInfo->setProteinRatio( queryProteinRatio );
}

vector< ProteinInfo * > ProteomeInfo::getProteinInfo( string sKeyword )
{
	if( sKeyword == "" )
		return vpProteinInfo;
	
	string::size_type positionLocus;
	string::size_type positionDescription;

	vector< ProteinInfo * > vpProteinInfoOutput;	
	int i;
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		positionLocus = vpProteinInfo[i]->getLocus().find( sKeyword );
		positionDescription = vpProteinInfo[i]->getDescription().find( sKeyword );

		if( positionLocus != string::npos || positionDescription != string::npos )
		{
			vpProteinInfoOutput.push_back( vpProteinInfo[i] );
		}
	}

	return vpProteinInfoOutput;
}

vector< ProteinInfo * > ProteomeInfo::getProteinInfo4Locus( string sLocus )
{
	vector< ProteinInfo * > vpProteinInfoOutput;	

	if( sLocus == "" )
		return vpProteinInfoOutput;

	int i;
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		if( sLocus == vpProteinInfo[i]->getLocus() )
		{
			vpProteinInfoOutput.push_back( vpProteinInfo[i] );
		}
	}

	return vpProteinInfoOutput;
}

void ProteomeInfo::getLocusList( vector< string > & vsLocusListRef )
{
	vsLocusListRef.clear();
	int i;
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		vsLocusListRef.push_back( vpProteinInfo[i]->getLocus() );
	}	

}

void ProteomeInfo::getLocusDescriptionList( vector< string > & vsLocusListRef , vector< string > & vsDescriptionRef)
{
	vsLocusListRef.clear();
	vsDescriptionRef.clear();
	int i;
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		vsLocusListRef.push_back( vpProteinInfo[i]->getLocus() );
		vsDescriptionRef.push_back( vpProteinInfo[i]->getDescription() );
	}	

}

void ProteomeInfo::sortProteinInfo( vector< ProteinInfo * > & vpProteinInfoInput, string sKey )
{	
	LessProteinInfo lessProteinInfoSort( sKey );
	sort( vpProteinInfoInput.begin(), vpProteinInfoInput.end(), lessProteinInfoSort );	
}

void ProteomeInfo::sortPeptideInfo( vector< PeptideInfo * > & vpPeptideInfoInput, string sKey )
{	
	LessPeptideInfo lessPeptideInfoSort( sKey );
	sort( vpPeptideInfoInput.begin(), vpPeptideInfoInput.end(), lessPeptideInfoSort );	
}

void ProteomeInfo::sortProteinInfoDescending( vector< ProteinInfo * > & vpProteinInfoInput, string sKey )
{	
	LessProteinInfo lessProteinInfoSort( sKey );
	sort( vpProteinInfoInput.begin(), vpProteinInfoInput.end(), lessProteinInfoSort );
	reverse( vpProteinInfoInput.begin(), vpProteinInfoInput.end() );
}

void ProteomeInfo::sortPeptideInfoDescending( vector< PeptideInfo * > & vpPeptideInfoInput, string sKey )
{	
	LessPeptideInfo lessPeptideInfoSort( sKey );
	sort( vpPeptideInfoInput.begin(), vpPeptideInfoInput.end(), lessPeptideInfoSort );
	reverse( vpPeptideInfoInput.begin(), vpPeptideInfoInput.end() );	
}	

bool ProteomeInfo::writeFileQPR()
{
	string sFilename = ProRataConfig::getWorkingDirectory() + "ProRata_Quantification.qpr.xml";
	FILE * pFileInitial;
	if( ( pFileInitial = fopen( sFilename.c_str(), "w" ) ) == NULL  ) 
	{
		cout << "ERROR: cannot write the result QPR file: " << sFilename << endl;
		return false;
	}

	// declaration
	ostringstream ossStream;
	ossStream << "<?xml version = \"1.0\" ?>" << "\n";
	ossStream << "<QUANTITATIVE_PROTEOMICS program = \"ProRata\"  version = \"" << ProRataConfig::getProRataVersion() << "\"> \n" ;
	fputs( ossStream.str().c_str(), pFileInitial );
	fclose( pFileInitial );	

	// write the config for peptide quantificationa and protein quantification to the QPR file
	if( !ProRataConfig::writeConfigXML( sFilename, 1, false, true, true ) )
	{
		cout << "ERROR: cannot write the CONFIG element to the result QPR file: " << sFilename << endl;
		return false;
	}
	
	// re-open the file in append mode
	FILE * pFILE;
	if( ( pFILE = fopen( sFilename.c_str(), "a" ) ) == NULL  ) 
	{
		cout << "ERROR: cannot write the result QPR file: " << sFilename << endl;
		return false;
	}
	
	int i;
	int j;
	int k;

	// sort the vpProteinInfo by their locus before writing to the file
	sortProteinInfo( vpProteinInfo, "locus" );
	
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		TiXmlElement * pElementProtein = new TiXmlElement( "PROTEIN" );
		
		// creat LOCUS
		TiXmlElement * pElementLocus = new TiXmlElement( "LOCUS" );
		TiXmlText * pTextLocus = new TiXmlText( vpProteinInfo[i]->getLocus().c_str() );
		pElementLocus->LinkEndChild( pTextLocus );
		pElementProtein->LinkEndChild( pElementLocus );

		// creat DESCRIPTION
		TiXmlElement * pElementDescription = new TiXmlElement( "DESCRIPTION" );
		TiXmlText * pTextDescription = new TiXmlText( vpProteinInfo[i]->getDescription().c_str() );
		pElementDescription->LinkEndChild( pTextDescription );
		pElementProtein->LinkEndChild( pElementDescription );

		// creat LOG2_RATIO
		TiXmlElement * pElementLog2Ratio = new TiXmlElement( "LOG2RATIO" );
		ossStream.str("");
		ossStream <<  vpProteinInfo[i]->getLog2Ratio(); 
		TiXmlText * pTextLog2Ratio = new TiXmlText( ossStream.str().c_str() );
		pElementLog2Ratio->LinkEndChild( pTextLog2Ratio );
		pElementProtein->LinkEndChild( pElementLog2Ratio );
		
		// creat CONFIDENCE_INTERVAL
		TiXmlElement * pElementCI = new TiXmlElement( "CONFIDENCE_INTERVAL" );

		TiXmlElement * pElementUpperLimitCI = new TiXmlElement( "UPPER_LIMIT" );
		ossStream.str("");
		ossStream <<  vpProteinInfo[i]->getUpperLimitCI(); 
		TiXmlText * pTextUpperLimitCI = new TiXmlText( ossStream.str().c_str() );
		pElementUpperLimitCI->LinkEndChild( pTextUpperLimitCI );
		pElementCI->LinkEndChild( pElementUpperLimitCI );
			
		TiXmlElement * pElementLowerLimitCI = new TiXmlElement( "LOWER_LIMIT" );
		ossStream.str("");
		ossStream <<  vpProteinInfo[i]->getLowerLimitCI(); 
		TiXmlText * pTextLowerLimitCI = new TiXmlText( ossStream.str().c_str() );
		pElementLowerLimitCI->LinkEndChild( pTextLowerLimitCI  );
		pElementCI->LinkEndChild(pElementLowerLimitCI );

		pElementProtein->LinkEndChild( pElementCI );

		vector< PeptideInfo * > vpTempPeptideInfo = vpProteinInfo[i]->getPeptideInfo();
		
		sortPeptideInfo( vpTempPeptideInfo, "sequence" );
		for( j = 0; j < vpTempPeptideInfo.size(); ++j )
		{
			TiXmlElement * pElementPeptide = new TiXmlElement( "PEPTIDE" );
			
			TiXmlElement * pElementChroFilename = new TiXmlElement( "XIC_FILE" );
			TiXmlText * pTextChroFilename = new TiXmlText( vpTempPeptideInfo[j]->getFilename().c_str() );
			pElementChroFilename->LinkEndChild( pTextChroFilename );
			pElementPeptide->LinkEndChild( pElementChroFilename );

			TiXmlElement * pElementIdentifier = new TiXmlElement( "XIC_IDENTIFIER" );
			ossStream.str("");
			ossStream <<  vpTempPeptideInfo[j]->getIdentifier(); 
			TiXmlText * pTextIdentifier = new TiXmlText( ossStream.str().c_str() );
			pElementIdentifier->LinkEndChild( pTextIdentifier );
			pElementPeptide->LinkEndChild( pElementIdentifier );

			TiXmlElement * pElementValidity = new TiXmlElement( "VALIDITY" );
			ossStream.str("");
			ossStream << boolalpha << vpTempPeptideInfo[j]->getValidity(); 
			TiXmlText * pTextValidity = new TiXmlText( ossStream.str().c_str() );
			pElementValidity->LinkEndChild( pTextValidity );
			pElementPeptide->LinkEndChild( pElementValidity );

			TiXmlElement * pElementPCALog2Ratio = new TiXmlElement( "LOG2RATIO" );
			ossStream.str("");
			ossStream << vpTempPeptideInfo[j]->getPCALog2Ratio(); 
			TiXmlText * pTextPCALog2Ratio = new TiXmlText( ossStream.str().c_str() );
			pElementPCALog2Ratio->LinkEndChild( pTextPCALog2Ratio );
			pElementPeptide->LinkEndChild( pElementPCALog2Ratio );

			TiXmlElement * pElementLog2SNR = new TiXmlElement( "LOG2SN" );
			ossStream.str("");
			ossStream << vpTempPeptideInfo[j]->getPCALog2SNR(); 
			TiXmlText * pTextPCALog2SNR = new TiXmlText( ossStream.str().c_str() );
			pElementLog2SNR->LinkEndChild( pTextPCALog2SNR );
			pElementPeptide->LinkEndChild( pElementLog2SNR );

			TiXmlElement * pElementSequence = new TiXmlElement( "SEQUENCE" );
			TiXmlText * pTextSequence = new TiXmlText( vpTempPeptideInfo[j]->getSequence().c_str() );
			pElementSequence->LinkEndChild( pTextSequence );
			pElementPeptide->LinkEndChild( pElementSequence );

			TiXmlElement * pElementChargeState = new TiXmlElement( "Z" );
			ossStream.str("");
			ossStream << vpTempPeptideInfo[j]->getChargeState(); 
			TiXmlText * pTextChargeState = new TiXmlText( ossStream.str().c_str() );			
			pElementChargeState->LinkEndChild( pTextChargeState );
			pElementPeptide->LinkEndChild( pElementChargeState );

			TiXmlElement * pElementScore = new TiXmlElement( "ID_SCORE" );
			ossStream.str("");
			ossStream << vpTempPeptideInfo[j]->getMaximumScore(); 
			TiXmlText * pTextScore = new TiXmlText( ossStream.str().c_str() );			
			pElementScore ->LinkEndChild( pTextScore );
			pElementPeptide->LinkEndChild( pElementScore );

			vector< string > vsLocus = vpTempPeptideInfo[j]->getLocus();
			for( k = 0; k < vsLocus.size(); ++k )
			{
				TiXmlElement * pElementLocus = new TiXmlElement( "LOCUS" );
				TiXmlText * pTextLocus = new TiXmlText( vsLocus[k].c_str() );
				pElementLocus->LinkEndChild( pTextLocus );
				pElementPeptide->LinkEndChild( pElementLocus );
			}

			pElementProtein->LinkEndChild( pElementPeptide );
		}	
		pElementProtein->Print( pFILE, 1 );
		fputs( "\n", pFILE );
		delete pElementProtein;
	}	
	
	ossStream.str("");
	ossStream << "</QUANTITATIVE_PROTEOMICS>" << "\n";
	fputs( ossStream.str().c_str(), pFILE );
	fclose( pFILE );			
	
	return true;

}

bool ProteomeInfo::readFileQPR( string sFilename )
{
	int i;
	int j;
	int k;

	for( i = 0; i < vpProteinInfo.size(); ++i )
		delete vpProteinInfo[i];
	vpProteinInfo.clear();

	for( j = 0; j < vpProteinInfo.size(); ++j )
		delete vpPeptideInfo[j];
	vpPeptideInfo.clear();

	// Creat a TinyXML document
	TiXmlDocument txdQPRFile;
	string sValue;
	istringstream issStream;
	
	// Try loading the file.
	if ( ! ( txdQPRFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading qpr.xml file" << sFilename << endl;
		return false;
	}

	TiXmlElement * pElementRoot = txdQPRFile.FirstChildElement( "QUANTITATIVE_PROTEOMICS" );
	
	vector<string> vsTagList;

	// start reading CONFIG
	vsTagList.clear();
	vsTagList.push_back( "CONFIG" );

	string sTemp;
	// extract elements inside PEPTIDE_QUANTIFICATION
	vsTagList.push_back( "PEPTIDE_QUANTIFICATION" );
	vsTagList.push_back( "PEAK_DETECTION" );
	vsTagList.push_back( "CHROMATOGRAM_SMOOTHING" );
	vsTagList.push_back( "WINDOW_SIZE" );
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	int iPeakDetectionWindowSize;
	issStream >> iPeakDetectionWindowSize;
	ProRataConfig::setPeakPickerWindowSize( iPeakDetectionWindowSize );
	
	vsTagList[4] = "ORDER";
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	int iPeakDetectionFOrder;
	issStream >> iPeakDetectionFOrder;
	ProRataConfig::setPeakPickerFOrder( iPeakDetectionFOrder );
	
	vsTagList[3] = "PEAK_SHIFT";
	vsTagList[4] = "LEFT";
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	float fLeftPeakShift;
	issStream >> fLeftPeakShift;
	ProRataConfig::setLeftPeakShift( fLeftPeakShift );

	vsTagList[4] = "RIGHT";
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	float fRightPeakShift;
	issStream >> fRightPeakShift;
	ProRataConfig::setRightPeakShift( fRightPeakShift );
	
	vsTagList.pop_back();
	vsTagList[2] = "LOG2_RATIO";
	vsTagList[3] = "MINIMUM";
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	double dPCAMinLog2Ratio;
	issStream >> dPCAMinLog2Ratio;
	ProRataConfig::setPCAMinLog2Ratio( dPCAMinLog2Ratio );

	vsTagList[3] = "MAXIMUM";
	sTemp = getValue( pElementRoot, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	double dPCAMaxLog2Ratio;
	issStream >> dPCAMaxLog2Ratio;
	ProRataConfig::setPCAMaxLog2Ratio( dPCAMaxLog2Ratio );

	vsTagList[2] = "ABUNDANCE_RATIO";
	vsTagList[3] = "NUMERATOR_ISOTOPOLOGUE";
	string sNumeratorIsotopologue = getValue( pElementRoot, vsTagList ) ;
	ProRataConfig::setNumeratorIsotopologue( sNumeratorIsotopologue );

	vsTagList[3] = "DENOMINATOR_ISOTOPOLOGUE";
	string sDenominatorIsotopologue = getValue( pElementRoot, vsTagList ) ;
	ProRataConfig::setDenominatorIsotopologue( sDenominatorIsotopologue );

	vsTagList.pop_back();
	vsTagList[2] = "LOG2_SNR_CUTOFF";
	sTemp =  getValue( pElementRoot, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	double dMinLog2SNR;
	issStream >> dMinLog2SNR;	
	ProRataConfig::setMinLog2SNR( dMinLog2SNR );
	
	vsTagList[2] = "REMOVE_AMBIGUOUS_PEPTIDES";
	sTemp = getValue( pElementRoot, vsTagList );
	bool bRemoveAmbiguousPeptide = true;
	if( sTemp == "true" || sTemp == "True" || sTemp == "TRUE" || sTemp == "T" || sTemp == "" )
		bRemoveAmbiguousPeptide = true;
	else
		bRemoveAmbiguousPeptide = false;
	ProRataConfig::setRemoveAmbiguousPeptides( bRemoveAmbiguousPeptide );
	
	// extract elements inside PROTEIN_QUANTIFICATION
	vsTagList.clear();
	vsTagList.push_back( "CONFIG" );
	vsTagList.push_back( "PROTEIN_QUANTIFICATION" );

	int iMinPeptideNumber = 2;
	issStream.clear();
	vsTagList.push_back( "MIN_PEPTIDE_NUMBER" );
	sTemp =  getValue( pElementRoot, vsTagList );
	if( sTemp != "" )
	{
		issStream.str( sTemp );
		issStream >> iMinPeptideNumber;
	}
	ProRataConfig::setMinPeptideNumber( iMinPeptideNumber );

	double dLog2RatioDiscretization = 0.1;
	issStream.clear();
	vsTagList[2] = "LOG2_RATIO_DISCRETIZATION";
	sTemp =  getValue( pElementRoot, vsTagList );
	if( sTemp != "" )
	{
		issStream.str( sTemp );
		issStream >> dLog2RatioDiscretization;
	}
	ProRataConfig::setLog2RatioDiscretization( dLog2RatioDiscretization );

	double dMaxCIwidth = 7.0;
	vsTagList[2] = "MAX_CI_WIDTH";
	sTemp =  getValue( pElementRoot, vsTagList );
	issStream.clear();
	if( sTemp != "" )
	{
		issStream.str( sTemp );
		issStream >> dMaxCIwidth;
	}
	ProRataConfig::setMaxCIwidth( dMaxCIwidth );
	
	double dMLEMinLog2Ratio = -7.0;
	vsTagList[2] = "LOG2_RATIO";
	vsTagList.push_back( "MINIMUM" );
	sTemp =  getValue( pElementRoot, vsTagList );
	issStream.clear();
	if( sTemp != "" )
	{
		issStream.str( sTemp );
		issStream >> dMLEMinLog2Ratio;	
	}
	ProRataConfig::setMLEMinLog2Ratio( dMLEMinLog2Ratio );

	double dMLEMaxLog2Ratio = 7.0;
	vsTagList[3] = "MAXIMUM";
	sTemp =  getValue( pElementRoot, vsTagList );
	issStream.clear();
	if( sTemp != "" )
	{
		issStream.str( sTemp );
		issStream >> dMLEMaxLog2Ratio;	
	}
	ProRataConfig::setMLEMaxLog2Ratio( dMLEMaxLog2Ratio );

	// start reading PROTEIN nodes
	vsTagList.clear();
	vsTagList.push_back( "PROTEIN" );

	vector< TiXmlElement * > vpElementVector;
	vector< TiXmlElement * > vpElementPeptideVector;
	vector< TiXmlElement * > vpElementLocusVector;
	
	vpElementVector = getElement( pElementRoot, vsTagList );


	for( i = 0; i < vpElementVector.size(); ++i )
	{
		ProteinInfo * pCurrentProteinInfo = new ProteinInfo;
		
		vsTagList.clear();
		vsTagList.push_back( "LOCUS" );
		sValue = getValue( vpElementVector[i], vsTagList );
		pCurrentProteinInfo->setLocus( sValue );
		
		vsTagList.clear();
		vsTagList.push_back( "DESCRIPTION" );
		sValue = getValue( vpElementVector[i], vsTagList );
		pCurrentProteinInfo->setDescription( sValue );
		
		vsTagList.clear();
		vsTagList.push_back( "LOG2RATIO" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		double dLog2Ratio;
		issStream >> dLog2Ratio;
		pCurrentProteinInfo->setLog2Ratio( dLog2Ratio );
		
		vsTagList.clear();
		vsTagList.push_back( "CONFIDENCE_INTERVAL" );
		vsTagList.push_back( "UPPER_LIMIT" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		double dUpperLimitCI;
		issStream >> dUpperLimitCI;
		pCurrentProteinInfo->setUpperLimitCI( dUpperLimitCI );

		vsTagList.pop_back();
		vsTagList.push_back( "LOWER_LIMIT" );
		sValue = getValue( vpElementVector[i], vsTagList );
		issStream.clear();
		issStream.str( sValue );
		double dLowerLimitCI;
		issStream >> dLowerLimitCI;
		pCurrentProteinInfo->setLowerLimitCI( dLowerLimitCI );	
		
		vsTagList.clear();
		vsTagList.push_back( "PEPTIDE" );		
		vpElementPeptideVector = getElement( vpElementVector[i], vsTagList );

		for( j = 0; j < vpElementPeptideVector.size(); ++j )
		{
			PeptideInfo * pCurrentPeptideInfo = new PeptideInfo;
			vsTagList.clear();
			vsTagList.push_back( "XIC_FILE" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			pCurrentPeptideInfo->setFilename( sValue );

			vsTagList.clear();
			vsTagList.push_back( "XIC_IDENTIFIER" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			issStream.clear();
			issStream.str( sValue );
			int iIdentifier;
			issStream >> iIdentifier;
			pCurrentPeptideInfo->setIdentifier( iIdentifier );

			vsTagList.clear();
			vsTagList.push_back( "VALIDITY" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			bool bValidity = false;
			if( sValue == "true" || sValue == "True" || sValue == "TRUE") 
				bValidity = true;
			pCurrentPeptideInfo->setValidity( bValidity );

			vsTagList.clear();
			vsTagList.push_back( "LOG2RATIO" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			issStream.clear();
			issStream.str( sValue );
			double dPCALog2Ratio;
			issStream >> dPCALog2Ratio;
			pCurrentPeptideInfo->setPCALog2Ratio( dPCALog2Ratio );			

			vsTagList.clear();
			vsTagList.push_back( "LOG2SN" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			issStream.clear();
			issStream.str( sValue );
			double dPCALog2SNR;
			issStream >> dPCALog2SNR;
			pCurrentPeptideInfo->setPCALog2SNR( dPCALog2SNR );

			vsTagList.clear();
			vsTagList.push_back( "SEQUENCE" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			pCurrentPeptideInfo->setSequence( sValue );

			vsTagList.clear();
			vsTagList.push_back( "Z" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			issStream.clear();
			issStream.str( sValue );
			int iChargeState;
			issStream >> iChargeState;
			pCurrentPeptideInfo->setChargeState( iChargeState );

			vsTagList.clear();
			vsTagList.push_back( "ID_SCORE" );
			sValue = getValue( vpElementPeptideVector[j], vsTagList );
			issStream.clear();
			issStream.str( sValue );
			float fScoreInput;
			issStream >> fScoreInput;
			pCurrentPeptideInfo->setMaximumScore( fScoreInput );
			
			vsTagList.clear();
			vsTagList.push_back( "LOCUS" );
			vpElementLocusVector = getElement( vpElementPeptideVector[j], vsTagList );
			vector< string > vsLocus;

			for( k = 0; k < vpElementLocusVector.size(); ++k )
			{
				vsTagList.clear();
				sValue = getValue( vpElementLocusVector[k], vsTagList );
				vsLocus.push_back( sValue );
			}
			pCurrentPeptideInfo->setLocus( vsLocus );		

			pCurrentProteinInfo->addPeptideInfo( pCurrentPeptideInfo );
			vpPeptideInfo.push_back( pCurrentPeptideInfo );
		}
	
		vpProteinInfo.push_back( pCurrentProteinInfo );
		
	}
	
	return true;
}

/*
 * An utility method to safely extract information from an XML tag.
 * the text inside an XML element is extracted and pasted together if separated
 * by comments or elements.
 * the element is reached by giving a vector of the element name in the order
 * of their hierarchy. Arbitrary number of level can be reached  
 */


string ProteomeInfo::getValue( TiXmlElement * pElement, const vector<string> &vsTagList )
{
	string sTemp = "";
	
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string
	if( !pElement )
		return sTemp;

	// creat a pointer to a node
	TiXmlNode * txnTemp = NULL;

	TiXmlText *txs;


	// check if the tree path is not empty
	if ( vsTagList.size() == 0 )
	{
		for( txnTemp = pElement->FirstChild(); txnTemp; txnTemp = txnTemp->NextSibling() )
		{
			// if this node is pointing to a node of type TEXT, which equals 4 in enum NodeType
			if( txnTemp->Type() == 4 )
			{
				// cast txnTemp to a text node
				txs = txnTemp->ToText();
				// get txnTemp's value and then append it to sTemp
				if( txs )
					sTemp.append( txs->Value() );
			}
		}

		return sTemp;

	}
	
	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	/*
	 * the lins is changed from this function's implementation for ProRataConfig
	 */
	// move the pointer to txddoc's first child with the specified tag name
	txnTemp = pElement->FirstChild( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! txnTemp )
	{
		cout << "ERROR: TAG\"" << (*itrTagListItr) << 
			"\" not found in the qpr file." << endl;
		return string("");
	}
	itrTagListItr++;
	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		txnTemp = txnTemp->FirstChild( (*itrTagListItr ).c_str() );

		if ( ! txnTemp )
		{
			cout << "ERROR: TAG\"" << (*itrTagListItr) << 
				"\" not found in the qpr file." << endl;
			return string("");
		}

	}

	// move the iterator back to point it to the last element name
	itrTagListItr--;

	/*
	 * inside the pointed element, there could be a mixture of
	 * text nodes, comment nodes and element nodes
	 * loop thru every nodes and for each text nodes, retrieve their text
	 * concatenate the text together and return them
	 */
	

	
	// point txnTemp to the child nodes and loop thru every child node
	for( txnTemp = txnTemp->FirstChild(); txnTemp; txnTemp = txnTemp->NextSibling() )
	{
		// if this node is pointing to a node of type TEXT, which equals 4 in enum NodeType
		if( txnTemp->Type() == 4 )
		{
			// cast txnTemp to a text node
			txs = txnTemp->ToText();
			// get txnTemp's value and then append it to sTemp
			if( txs )
				sTemp.append( txs->Value() );
		}
	}

	return sTemp;
	


}

vector< TiXmlElement * > ProteomeInfo::getElement(  TiXmlElement * pElement, const vector<string> &vsTagList )
{
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string

	vector< TiXmlElement * > vpElementVector;

	if( !pElement )
		return vpElementVector;

	// check if the tree path is not empty
	if ( vsTagList.size() < 1 )
		return vpElementVector;
	
	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	/*
	 * the lins is changed from this function's implementation for ProRataConfig
	 */
	// move the pointer to txddoc's first child with the specified tag name
	pElement = pElement->FirstChildElement( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! pElement )
	{
		cout << "ERROR: TAG\"" << (*itrTagListItr) << 
			"\" not found in the qpr file." << endl;
		return vpElementVector;
	}

	itrTagListItr++;

	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		pElement = pElement->FirstChildElement( (*itrTagListItr).c_str() );

		if ( ! pElement )
		{
			cout << "ERROR: TAG\"" << (*itrTagListItr) << 
				"\" not found in the qpr file." << endl;
			return vpElementVector;
		}

	}

	itrTagListItr--;

	TiXmlElement * pElementCurrent = pElement;
	for( pElement = pElementCurrent; pElement; pElement = pElement->NextSiblingElement( (*itrTagListItr).c_str() ) )
	{
		vpElementVector.push_back( pElement );
	}

	return vpElementVector;
}

bool ProteomeInfo::writeFileTAB()
{
	unsigned int i;
	// write the flat-file table for the protein quantification results
	string sCompleteFilename = ProRataConfig::getWorkingDirectory() + "ProRata_Quantification_Protein.txt";
	ofstream fStreamTro( sCompleteFilename.c_str() );
	if( !fStreamTro )
	{
		cout << "ERROR: cannot write ProRata_Quantification_Protein.txt" << endl;
		return false;
	}

	fStreamTro << "Protein Quantification Results by ProRata " << ProRataConfig::getProRataVersion() <<  endl;
	fStreamTro << endl;
	fStreamTro << "Configurations" << endl;
	fStreamTro << "Ratio of samples = " << ProRataConfig::getNumeratorIsotopologue() << " : " <<  ProRataConfig::getDenominatorIsotopologue() << endl;
	fStreamTro << "Minimum number of quantified peptides = " << ProRataConfig::getMinPeptideNumber() << endl;
	fStreamTro << "Maximum confidence interval width cutoff = " << ProRataConfig::getMaxCIwidth() << endl;	
	fStreamTro << "Range of the log2 protein abundance ratio = [" << ProRataConfig::getMLEMinLog2Ratio() << ", " <<  ProRataConfig::getMLEMaxLog2Ratio() << "]" << endl;
	fStreamTro << "Discretization precision of the log2 ratio = " << ProRataConfig::getLog2RatioDiscretization() << endl;
	fStreamTro << "Ln-likelihood cutoff offset = " << ProRataConfig::getLnLikelihoodCutoffOffset() << endl;
	fStreamTro << endl;
	fStreamTro << "Tab-delimited Table" << endl;
	fStreamTro << "locus" << '\t' 
		<< "log2ratio" << '\t'
		<< "lowerCI" << '\t'
		<< "upperCI" << '\t'
		<< "quantified_peptides" << '\t'
		<< "description" << endl;
	sortProteinInfo( vpProteinInfo, "locus" );
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		fStreamTro << vpProteinInfo[i]->getLocus() << '\t' 
			<< vpProteinInfo[i]->getLog2Ratio() << '\t'
			<< vpProteinInfo[i]->getLowerLimitCI() << '\t'
			<< vpProteinInfo[i]->getUpperLimitCI() << '\t'
			<< vpProteinInfo[i]->getQuantifiedPeptides() << '\t'
			<< vpProteinInfo[i]->getDescription() << endl;
	}
	
	fStreamTro.close();

	sCompleteFilename = ProRataConfig::getWorkingDirectory() +  "ProRata_Quantification_Peptide.txt";
	// write the flat-file table for the peptide quantification results
	ofstream fStreamTep( sCompleteFilename.c_str()  );
	if( !fStreamTep )
	{
		cout << "ERROR: cannot write ProRata_Quantification_Peptide.txt" << endl;
		return false;
	}

	fStreamTep << "Peptide Quantification Results by ProRata " << ProRataConfig::getProRataVersion() << endl;
	fStreamTep << endl;
	fStreamTep << "Configurations" << endl;
	fStreamTep << "Ratio of samples = " << ProRataConfig::getNumeratorIsotopologue() << " : " <<  ProRataConfig::getDenominatorIsotopologue() << endl;
	fStreamTep << "Range of the log2 peptide abundance ratio = [" << ProRataConfig::getPCAMinLog2Ratio() << ", " <<  ProRataConfig::getPCAMaxLog2Ratio() << "]" << endl;
	fStreamTep << endl;
	fStreamTep << "Tab-delimited Table" << endl;
	fStreamTep << "XIC_filename" << '\t' 
		<< "identifier" << '\t'
		<< "log2Ratio" << '\t'
		<< "log2SNR" << '\t'
		<< "validity" << '\t'
		<< "locus" << '\t'
		<< "sequence" << endl;
	for( i = 0; i < vpPeptideInfo.size(); ++i )
	{
		fStreamTep << vpPeptideInfo[i]->getFilename() << '\t' 
			<< vpPeptideInfo[i]->getIdentifier() << '\t'
			<< vpPeptideInfo[i]->getPCALog2Ratio() << '\t'
			<< vpPeptideInfo[i]->getPCALog2SNR() << '\t'
			<< boolalpha << vpPeptideInfo[i]->getValidity() << '\t';
		
		vector< string > vsLocus = vpPeptideInfo[i]->getLocus();
		for(int k = 0; k <  vsLocus.size(); k++ )
		{
			fStreamTep << vsLocus[k];
			if( k != (vsLocus.size() - 1) )
			{
				fStreamTep << ",";
			}
				
		}
		fStreamTep << '\t' << vpPeptideInfo[i]->getSequence() << endl;
	}

	fStreamTep.close();
	return true;
}


bool ProteomeInfo::writeFileLabelFree()
{
	unsigned int i;
	// write the flat-file table for the protein quantification results
	string sCompleteFilename = ProRataConfig::getWorkingDirectory() + "ProRata_LabelFree_Quantification_Protein.txt";
	ofstream fStreamTro( sCompleteFilename.c_str() );
	if( !fStreamTro )
	{
		cout << "ERROR: cannot write ProRata_LabelFree_Quantification_Protein.txt " << endl;
		return false;
	}

	fStreamTro << "Protein Quantification Results by ProRata " << ProRataConfig::getProRataVersion() <<  endl;
	fStreamTro << "locus" << '\t' 
		<< "totalPeakHeight" << '\t'
		<< "totalPeakArea" << '\t'
		<< "MS2Count" << '\t'
		<< "quantified_peptides" << '\t'
		<< "description" << endl;
	sortProteinInfo( vpProteinInfo, "locus" );
	for( i = 0; i < vpProteinInfo.size(); ++i )
	{
		fStreamTro << vpProteinInfo[i]->getLocus() << '\t' 
			<< vpProteinInfo[i]->getTotalPeakHeight() << '\t'
			<< vpProteinInfo[i]->getTotalPeakArea() << '\t'
			<< vpProteinInfo[i]->getMS2SpectralCounts() << '\t'
			<< vpProteinInfo[i]->getQuantifiedPeptides() << '\t'
			<< vpProteinInfo[i]->getDescription() << endl;
	}
	
	fStreamTro.close();

	sCompleteFilename = ProRataConfig::getWorkingDirectory() +  "ProRata_LabelFree_Quantification_Peptide.txt";
	// write the flat-file table for the peptide quantification results
	ofstream fStreamTep( sCompleteFilename.c_str()  );
	if( !fStreamTep )
	{
		cout << "ERROR: cannot write ProRata_LabelFree_Quantification_Peptide.txt" << endl;
		return false;
	}

	fStreamTep << "Peptide Quantification Results by ProRata " << ProRataConfig::getProRataVersion() << endl;
	fStreamTep << "XIC_filename" << '\t' 
		<< "identifier" << '\t'
		<< "peakHeight" << '\t'
		<< "peakArea" << '\t'
		<< "peakSNR" << '\t'
		<< "peakWidth" << '\t'
		<< "peakStartRT" << '\t'
		<< "peakEndRT" << '\t'
		<< "MS2Count" << '\t'
		<< "MS2RT" << '\t'
		<< "validity" << '\t'
		<< "locus" << '\t'
		<< "sequence" << '\t'
		<< "charge" << '\t'
		<< "IDs" << endl;
	for( i = 0; i < vpPeptideInfo.size(); ++i )
	{
		// build MS2Time
		stringstream ssMS2Time;
		vector<float> vfMS2Time = vpPeptideInfo[i]->getMS2Time();	
		for(int j = 0; j < vfMS2Time.size() - 1; ++j )
		{
			ssMS2Time << vfMS2Time[j] << ",";
		}
		ssMS2Time << vfMS2Time.back();

		fStreamTep << vpPeptideInfo[i]->getFilename() << '\t' 
			<< vpPeptideInfo[i]->getIdentifier() << '\t'
			<< vpPeptideInfo[i]->getPeakHeight() << '\t'
			<< vpPeptideInfo[i]->getPeakArea() << '\t'
			<< vpPeptideInfo[i]->getPeakSNR() << '\t'
			<< vpPeptideInfo[i]->getPeakTimeWidth() << '\t'
			<< vpPeptideInfo[i]->getLeftValleyTime() << '\t'
			<< vpPeptideInfo[i]->getRightValleyTime() << '\t'			
			<< vpPeptideInfo[i]->getMS2Count() << '\t'
			<< ssMS2Time.str() << '\t'
			<< boolalpha << vpPeptideInfo[i]->getValidity() << '\t';
		
		// write Locus
		vector< string > vsLocus = vpPeptideInfo[i]->getLocus();
		for(int k = 0; k <  vsLocus.size(); k++ )
		{
			fStreamTep << vsLocus[k];
			if( k != (vsLocus.size() - 1) )
			{
				fStreamTep << ",";
			}
				
		}
		fStreamTep << '\t' << vpPeptideInfo[i]->getSequence() ;
		fStreamTep << '\t' << vpPeptideInfo[i]->getChargeState() << '\t';

		vector<string> vsIDFilename = vpPeptideInfo[i]->getAllIDfilename();	
		for(int j = 0; j < vsIDFilename.size() - 1; ++j )
		{
			fStreamTep << vsIDFilename[j] << ",";
		}
		fStreamTep << vsIDFilename.back() << endl;
	}

	fStreamTep.close();
	return true;
}













