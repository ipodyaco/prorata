#include "isotopologue.h"

IsotopeDistribution::IsotopeDistribution()
{

}


IsotopeDistribution::IsotopeDistribution( vector< double > vItsMass,  vector< double > vItsProb )
{
	vMass = vItsMass;
	vProb = vItsProb;
}

IsotopeDistribution::~IsotopeDistribution()
{
	// destructor	
}

void IsotopeDistribution::print()
{
	cout << "Mass " << '\t' << "Inten" << endl;
	for( unsigned int i = 0; i < vMass.size(); i++)
	{
		cout << setprecision(8) << vMass[i] << '\t' << vProb[i] << endl;
	}
}

double IsotopeDistribution::getMostAbundantMass()
{
	double dMaxProb = 0;
	double dMass = 0;	
	for( unsigned int i = 0; i < vMass.size(); ++i )
	{
		if( dMaxProb < vProb[i] )
		{
			dMaxProb = vProb[i];
			dMass = vMass[i];
		}
	}
	return dMass;
}


double IsotopeDistribution::getAverageMass()
{
	double dSumProb = 0;
	double dSumMass = 0;	
	for( unsigned int i = 0; i < vMass.size(); ++i )
	{
		dSumProb = dSumProb + vProb[i];
		dSumMass = dSumMass + vProb[i] * vMass[i];
	}

	if(dSumProb <= 0)
		return 1.0;
	
	return ( dSumMass / dSumProb ) ;
}


Isotopologue::Isotopologue() : MassPrecision(0.1), ProbabilityCutoff(0.000000001), AtomName("CHONPS"), AtomNumber( 12 )
{
}

Isotopologue::Isotopologue( string sName,  string sTable ) : MassPrecision(0.1), ProbabilityCutoff(0.000000001), AtomName("CHONPS"), AtomNumber( 12 )
{

	setupIsotopologue( sName, sTable );
}

Isotopologue::~Isotopologue()
{
	// destructor
}

bool Isotopologue::setupIsotopologue( const string & sName,  const string & sTable )
{
	sIsotoplogueName = sName;
	
	istringstream issStream( sTable );
	string sResidue;
	vector< int > viAtomVector;
	int iNumber;
	int i;
	
	// parse out the RESIDUE_ATOMIC_COMPOSITION table
	while( !( issStream.eof() ) )
	{
		// each row is expected to start with the residue name, following by 12 numbers for the natuual and enriched CHONPS 
		issStream >> sResidue;
		if( sResidue == "" )
			continue;
		viAtomVector.clear();
		viAtomVector.reserve( AtomNumber );
		for( i = 0; i < AtomNumber; ++i )
		{
			if( issStream.eof() )
			{
				// this row doesn't have 12 fields
				cout << "ERROR:  the RESIDUE_ATOMIC_COMPOSITION table in ProRataConfig is not correct!" << endl;
				return false;
			}
			issStream >> iNumber;
			viAtomVector.push_back( iNumber );

		}
		// add this row into the mResidueAtomicComposition table
		mResidueAtomicComposition[ sResidue ] = viAtomVector;
		sResidue = "";
	}

	// push 12 empty IsotopeDistributions into vAtomIsotopicDistribution
	vAtomIsotopicDistribution.reserve( AtomNumber );
	for( i = 0; i < ( AtomNumber ); ++i )
	{
		IsotopeDistribution TempDistribution;
		vAtomIsotopicDistribution.push_back( TempDistribution );
	}
	
	// variables to be passed as reference to ProRataConfig::getAtomIsotopicComposition
	// to receive its return value
	vector< double > vdMassTemp;
	vector< double > vdNaturalCompositionTemp;
	vector< double > vdEnrichedCompositionTemp;
		
	// the isotopic distribution is pushed into vAtomIsotopicDistribution in the order of
	// natural CHONPS and then enriched CHONPS
	for( i = 0; i < (int)AtomName.size(); ++i )
	{
		if( ! ProRataConfig::getAtomIsotopicComposition( 
					AtomName[i],
					vdMassTemp,
					vdNaturalCompositionTemp,
					vdEnrichedCompositionTemp ) )
		{
			cout << "ERROR: cannot retrieve isotopic composition for atom " << AtomName[i] << " from ProRataConfig" << endl;
			return false;
		}
		vAtomIsotopicDistribution[i].vMass = vdMassTemp;
		vAtomIsotopicDistribution[i].vProb = vdNaturalCompositionTemp;
		vAtomIsotopicDistribution[ ( i+AtomName.size() ) ].vMass = vdMassTemp;
		vAtomIsotopicDistribution[ ( i+AtomName.size() ) ].vProb = vdEnrichedCompositionTemp;
	}


	// calculate Isotopic distribution for all residues
	map< string, vector< int > >::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for( ResidueIter = mResidueAtomicComposition.begin(); ResidueIter != mResidueAtomicComposition.end(); ResidueIter++ )
	{
		if( !computeIsotopicDistribution( ResidueIter->second, tempIsotopeDistribution ) )
		{
			cout << "ERROR: cannot calculate the isotopic distribution for residue " << ResidueIter->first << endl;
			return false;
		}
		vResidueIsotopicDistribution[ ResidueIter->first ] = tempIsotopeDistribution;
	}
	




	return true;

}

bool Isotopologue::computeMZwindows( string sSequence, int iChargeState, MZwindows & myMZwindows )
{
	
	formatSequence( sSequence ); 

	IsotopeDistribution distributionPeptide;
	if( !computeIsotopicDistribution( sSequence, distributionPeptide ) )
	{
		cout << "WARNNING: cannot compute isotopic distribution for peptide " << sSequence << endl;
		return false;
	}

	unsigned int i = 0;

	// get the masses that passes the isotopic envelop cutoff
	vector< double > vMassFiltered;
	double dMaxProb = *max_element( distributionPeptide.vProb.begin(), distributionPeptide.vProb.end() );
	double dCutoffProb = dMaxProb * ProRataConfig::getIsotopicEnvelopCutoff();
	for( i = 0; i < distributionPeptide.vProb.size(); i++ )
	{
		if( distributionPeptide.vProb[i] > dCutoffProb )
			vMassFiltered.push_back( distributionPeptide.vMass[i] );
	}

	if( iChargeState <= 0 )
	{
		cout << "ERROR: the charge state cannot be zero or negative number! " << endl;
		return false;
	}

	// calculate the mass to charge ratios
	vector< double > vMassToCharge;
	for( i = 0; i < vMassFiltered.size(); i++ )
		vMassToCharge.push_back( (vMassFiltered[i]+iChargeState)/((double)iChargeState) );
	
	// calculate the m/z ranges for each m/z value
	vector< double > vPeakUpperMZ;
	vector< double > vPeakLowerMZ;
	for( i = 0; i < vMassToCharge.size(); i++ )
	{
		vPeakLowerMZ.push_back( (vMassToCharge[i] - ProRataConfig::getMinusMZerror() ) );
		if( vPeakLowerMZ[i] < 0 )
			vPeakLowerMZ[i] = 0;
		vPeakUpperMZ.push_back( (vMassToCharge[i] + ProRataConfig::getPlusMZerror() ) );
	}

	// merge overlaping m/z ranges in vPeakUpperMZ and vPeakLowerMZ
	unsigned int j = 0;
	for( i = 0; i < vPeakUpperMZ.size(); i++ )
	{
		// determine if this m/z range has been merged previously
		if( (vPeakUpperMZ[i] > 0) && (vPeakLowerMZ[i] > 0) )
		{
			// compare with other m/z ranges
			for( j = (i+1); j < vPeakUpperMZ.size(); j++ )
			{
				if( (vPeakUpperMZ[j] > 0) && (vPeakLowerMZ[j] > 0) )
				{
					// determined if these two m/z ranges overlap
					if( ( -0.05 <= (vPeakUpperMZ[j] - vPeakLowerMZ[i]) ) && ( -0.05 <= (vPeakUpperMZ[i] - vPeakLowerMZ[j]) ) )
					{
						// if so, set the joint m/z range to m/z range [i] and flag m/z range [j] as having been merged 
						vPeakUpperMZ[i] = maximum( vPeakUpperMZ[i], vPeakUpperMZ[j] );
						vPeakLowerMZ[i] = minimum( vPeakLowerMZ[i], vPeakLowerMZ[j] );
						vPeakUpperMZ[j] = (-1.0);
						vPeakLowerMZ[j] = (-1.0);
					}
				}

			}	
		}
	}

	// save the merged m/z ranges to myMZwindows to be returned
	myMZwindows.vfUpperMZ.clear();
	myMZwindows.vfLowerMZ.clear();	
	for( i = 0; i < vPeakUpperMZ.size(); i++ )
	{
		if( (vPeakUpperMZ[i] > 0) && (vPeakLowerMZ[i] > 0) )
		{
			myMZwindows.vfUpperMZ.push_back( (float)vPeakUpperMZ[i] );
			myMZwindows.vfLowerMZ.push_back( (float)vPeakLowerMZ[i] );
		}
	}

	return true;

}

double Isotopologue::computeMostAbundantMass( string sSequence )
{
	IsotopeDistribution tempIsotopeDistribution;
	if( !computeIsotopicDistribution( sSequence, tempIsotopeDistribution ) )
		return 0;
	else
		return tempIsotopeDistribution.getMostAbundantMass();
}

double Isotopologue::computeAverageMass( string sSequence )
{
	IsotopeDistribution tempIsotopeDistribution;
	if( !computeIsotopicDistribution( sSequence, tempIsotopeDistribution ) )
		return 0;
	else
		return tempIsotopeDistribution.getAverageMass();
}

bool Isotopologue::computeIsotopicDistribution( string sSequence , IsotopeDistribution & myIsotopeDistribution )
{
	formatSequence( sSequence );

	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;
	map< string, IsotopeDistribution >::iterator ResidueIter;
	
	ResidueIter = vResidueIsotopicDistribution.find("NTerm");
	if(ResidueIter != vResidueIsotopicDistribution.end() )
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cout << "ERROR: can't find the N-terminus" << endl;
		return false;
	}
 
	ResidueIter = vResidueIsotopicDistribution.find("CTerm");
	if(ResidueIter != vResidueIsotopicDistribution.end() )
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = sum( currentDistribution, sumDistribution );
	}
	else
	{
		cout << "ERROR: can't find the C-terminus" << endl;
		return false;
	}

	// add up all residues's isotopic distribution
	for( unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr( j, 1 );
		ResidueIter = vResidueIsotopicDistribution.find( currentResidue );
		if(ResidueIter != vResidueIsotopicDistribution.end() )
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum( currentDistribution, sumDistribution );
		}
		else
		{
			cout << "ERROR: can't find the residue" << currentResidue << endl;
			return false;
		}
	}

	myIsotopeDistribution = sumDistribution;

	return true;

}

bool Isotopologue::computeProductIonMass( string sSequence, vector< double > & vdYion, vector< double > & vdBion )
{
	formatSequence( sSequence );
	
//	cout << "sSequence = " << sSequence << endl;
	
	vdYion.clear();
	vdBion.clear();

	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentDistribution;
	map< string, IsotopeDistribution >::iterator ResidueIter;	
	
	// compute B-ion series
	ResidueIter = vResidueIsotopicDistribution.find("NTerm");
	if(ResidueIter != vResidueIsotopicDistribution.end() )
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cout << "ERROR: can't find the N-terminus" << endl;
		return false;
	}
	
	for(unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr( j, 1 );
		ResidueIter = vResidueIsotopicDistribution.find( currentResidue );
		if(ResidueIter != vResidueIsotopicDistribution.end() )
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum( currentDistribution, sumDistribution );
		}
		else
		{
			cout << "ERROR: can't find the residue" << currentResidue << endl;
			return false;
		}
		// if the current residue is an alphabet then it is a amino acid residue ( not a PTM )
		// then save it to vdBion
		if( isalpha( currentResidue[0] ) )
			vdBion.push_back( sumDistribution.getAverageMass() );
	}
	
	// compute Y-ion series
	ResidueIter = vResidueIsotopicDistribution.find("CTerm");
	if(ResidueIter != vResidueIsotopicDistribution.end() )
	{
		currentDistribution = ResidueIter->second;
		sumDistribution = currentDistribution;
	}
	else
	{
		cout << "ERROR: can't find the N-terminus" << endl;
		return false;
	}
	
	for( int j = (sSequence.length() - 1) ; j > -1 ; j--)
	{
		string currentResidue = sSequence.substr( j, 1 );
		ResidueIter = vResidueIsotopicDistribution.find( currentResidue );
		if(ResidueIter != vResidueIsotopicDistribution.end() )
		{
			currentDistribution = ResidueIter->second;
			sumDistribution = sum( currentDistribution, sumDistribution );
		}
		else
		{
			cout << "ERROR: can't find the residue" << currentResidue << endl;
			return false;
		}
		// if the current residue is an alphabet then it is a amino acid residue ( not a PTM )
		// then save it to vdBion
		if( isalpha( currentResidue[0] ) )
			vdYion.push_back( ( sumDistribution.getAverageMass() + 2.0 ) );
		
	}
	return true;
	
}

bool Isotopologue::computeIsotopicDistribution( vector< int > AtomicComposition, IsotopeDistribution & myIsotopeDistribution )
{
	IsotopeDistribution sumDistribution;
	IsotopeDistribution currentAtomDistribution;
	currentAtomDistribution = multiply( vAtomIsotopicDistribution[0], AtomicComposition[0] );
	sumDistribution = currentAtomDistribution;

	for( int i = 1; i < AtomNumber; i++ )
	{
		currentAtomDistribution = multiply( vAtomIsotopicDistribution[i], AtomicComposition[i] );
		sumDistribution = sum( currentAtomDistribution, sumDistribution );
	}
	myIsotopeDistribution = sumDistribution;



	return true;
}

bool Isotopologue::computeAtomicComposition( string sSequence, vector< int >  & myAtomicComposition )
{
	formatSequence( sSequence );
	vector< int > AtomicComposition;
	vector< int > CurrentComposition;
	int i;
	map< string, vector< int > >::iterator ResidueIter;
	
	for( i = 0; i < AtomNumber; i++)
		AtomicComposition.push_back( 0 );
	
	ResidueIter = mResidueAtomicComposition.find("NTerm");
	if(ResidueIter != mResidueAtomicComposition.end() )
	{
		CurrentComposition = ResidueIter->second;
		for( i = 0; i < AtomNumber; i++)
			AtomicComposition[i] = CurrentComposition[i];
	}
	else
	{
		cout << "ERROR: can't find the atomic composition for the N-terminus" << endl;
		return false;
	}

	ResidueIter = mResidueAtomicComposition.find("CTerm");
	if(ResidueIter != mResidueAtomicComposition.end() )
	{
		CurrentComposition = ResidueIter->second;
		for( i = 0; i < AtomNumber; i++)
			AtomicComposition[i] += CurrentComposition[i];
	}
	else
	{
		cout << "ERROR: can't find the atomic composition for the C-terminus" << endl;
		return false;
	}


	for(unsigned int j = 0; j < sSequence.length(); j++)
	{
		string currentResidue = sSequence.substr( j, 1 );
		ResidueIter = mResidueAtomicComposition.find( currentResidue );
		if(ResidueIter != mResidueAtomicComposition.end() )
		{
			CurrentComposition = ResidueIter->second;
			for( i = 0; i < AtomNumber; i++)
				AtomicComposition[i] += CurrentComposition[i];
		}
		else
		{
			cout << "ERROR: can't find the atomic composition for residue/PTM: " << currentResidue << endl;	
			return false;
		}
	}


	myAtomicComposition = AtomicComposition;
	return true;

}

IsotopeDistribution Isotopologue::sum( IsotopeDistribution distribution0, IsotopeDistribution distribution1)
{


	IsotopeDistribution sumDistribution;
	double currentMass;
	double currentProb;
	bool bIsMerged;
	int iSizeDistribution0 = distribution0.vMass.size();
	int iSizeDistribution1 = distribution1.vMass.size();
	int iSizeSumDistribution;
	int i;
	int j;
	int k;
	for( i = 0; i < iSizeDistribution0; ++i )
	{
		for( j = 0; j < iSizeDistribution1; ++j )
		{
			// combine one isotopologue from distribution0 with one isotopologue distribution1
			currentMass = (distribution0.vMass[i] + distribution1.vMass[j]);
			currentProb = (distribution0.vProb[i] * distribution1.vProb[j]);
			if ( ( currentProb > ProbabilityCutoff )  ) 
			{
				iSizeSumDistribution = sumDistribution.vMass.size();
				// push back the first peak
				if( iSizeSumDistribution == 0 )
				{
					sumDistribution.vMass.push_back( currentMass );
					sumDistribution.vProb.push_back( currentProb );
				}
				else
				{
					bIsMerged = false;			
					// check if the combined isotopologue can be merged with the existing isotopologues 
					for ( k = 0; k < iSizeSumDistribution; ++k )
					{
						if( fabs( currentMass - sumDistribution.vMass[k] ) < MassPrecision )
						{
							// average Mass
							sumDistribution.vMass[k] = ( currentMass*currentProb + sumDistribution.vMass[k]*sumDistribution.vProb[k] ) / ( currentProb + sumDistribution.vProb[k] );
							// sum Prob
							sumDistribution.vProb[k] = currentProb + sumDistribution.vProb[k];
							bIsMerged = true;
							break;
						}
					}
					if( !bIsMerged )
					{
						sumDistribution.vMass.push_back( currentMass );
						sumDistribution.vProb.push_back( currentProb );
					}
				}
			}
		}
	}


	// normalize the probability space to 1
	double sumProb = 0;
	iSizeSumDistribution = sumDistribution.vMass.size();
	for( i = 0; i < iSizeSumDistribution; ++i )
		sumProb += sumDistribution.vProb[i];
	
	if( sumProb <= 0 )
		return sumDistribution;
	
	for( i = 0; i < iSizeSumDistribution; ++i )
		sumDistribution.vProb[i] = sumDistribution.vProb[i]/sumProb;
	
	return sumDistribution;

}


IsotopeDistribution Isotopologue::multiply( IsotopeDistribution distribution0, int count )
{
	if( count == 1 )
		return distribution0;
	IsotopeDistribution productDistribution;
	productDistribution.vMass.push_back( 0.0 );
	productDistribution.vProb.push_back( 1.0 );
	for( int i = 0; i < ( abs(count) ); ++i )
	{
		productDistribution = sum( productDistribution, distribution0 );
	}
	return productDistribution;

}

void Isotopologue::formatSequence( string & sSequence )
{
	int iLength = sSequence.length();
	// check if this sequence is of the format X.XXXX.X
	// its minimum length is 5
	if(iLength > 4 )
	{
		// its second residues from both the left side and the right side should be '.'
		if( sSequence[1] == '.' && sSequence[ ( iLength - 2 ) ] == '.')
		{
			// if so, extract the sequence between '.'
			sSequence = sSequence.substr( 2, (iLength - 4) );
		}
	}
}
	
