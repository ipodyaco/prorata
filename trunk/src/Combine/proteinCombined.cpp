
#include "proteinCombined.h"

ProteinCombined::ProteinCombined()
{
	sLocus = "";
	sDescription = "";
	bValidity = false;
}

ProteinCombined::~ProteinCombined()
{
	// destructor
}

ProteinCombined::ProteinCombined( const ProteinCombined & templateProteinCombined )
{
	vDirectComparison = templateProteinCombined.vDirectComparison;
	vIndirectComparison = templateProteinCombined.vIndirectComparison;
}

bool ProteinCombined::runQuantification( string sInputLocus )
{
	return runQuantification( sInputLocus, "" );
}


bool ProteinCombined::runQuantification( string sInputLocus, string sInputDescription )
{
	sLocus = sInputLocus;
	sDescription = sInputDescription;

	int i;
	for( i = 0; i < vDirectComparison.size(); ++i )
	{
		if( !vDirectComparison[i].runQuantification( sLocus ) )
		{
			cout << "WARNNING: problem in computing direct comparison for locus " << sLocus << endl;
		}
	}

	for( i = 0; i < vIndirectComparison.size(); ++i )
	{
		ProteinDirectComparison tempDirectComparison0;
		ProteinDirectComparison tempDirectComparison1;
		getDirectComparison( vIndirectComparison[i].getDirectComparisonName0(), tempDirectComparison0 ); 
		getDirectComparison( vIndirectComparison[i].getDirectComparisonName1(),  tempDirectComparison1 );

		if( !vIndirectComparison[i].runQuantification( sLocus, tempDirectComparison0, tempDirectComparison1 ) )
		{
			cout << "WARNNING: problem in computing indirect comparison for locus " << sLocus << endl;
		}
	}
	
	computeValidity();

	for( i = 0; i < vDirectComparison.size(); ++i )
	{
		vDirectComparison[i].eraseProfileLikelihoodCurve();
	}

	for( i = 0; i < vIndirectComparison.size(); ++i )
	{
		vIndirectComparison[i].eraseProfileLikelihoodCurve();
	}

	return true;
}

void ProteinCombined::eraseProfileLikelihoodCurve()
{
	int i;
	for( i = 0; i < vDirectComparison.size(); ++i )
	{
		vDirectComparison[i].eraseProfileLikelihoodCurve();
	}
	for( i = 0; i < vIndirectComparison.size(); ++i )
	{
		vIndirectComparison[i].eraseProfileLikelihoodCurve();
	}

}

string ProteinCombined::getLocus()
{
	return sLocus;
}

string ProteinCombined::getDescription()
{
	return sDescription;
}

bool ProteinCombined::getValidity()
{
	return bValidity;
}

bool ProteinCombined::addDirectComparison(ProteinDirectComparison tempProteinDirectComparison)
{
	// check  NUMERATOR and DENOMINATOR
	if( tempProteinDirectComparison.getNumerator() == tempProteinDirectComparison.getDenominator() )
	{
		cout << "ERROR: A comparison MUST have different names for its NUMERATOR and DENOMINATOR. Rename " 
			<< tempProteinDirectComparison.getNumerator() << endl;
		return false;
	}
		
	// check if this direct comparison has a name the same as the existing ones; 
	string sName = tempProteinDirectComparison.getName();
	if( IsThereDirectComparison( sName ) )
	{
		cout << "ERROR: A comparison MUST have a unique name! this name is used more than once, " << sName << endl;
		return false;
	}
	
	vDirectComparison.push_back( tempProteinDirectComparison );
	return true;
}

bool ProteinCombined::addIndirectComparison(ProteinIndirectComparison tempProteinIndirectComparison)
{
	// check NUMERATOR and DENOMINATOR
	if( tempProteinIndirectComparison.getNumerator() == tempProteinIndirectComparison.getDenominator() )
	{
		cout << "ERROR: A comparison MUST have different names for its NUMERATOR and DENOMINATOR. Rename " 
			<< tempProteinIndirectComparison.getNumerator() << endl;
		return false;
	}
		
	// check if this indirect comparison has a name the same as the existing ones; 
	string sName = tempProteinIndirectComparison.getName();
	if( IsThereIndirectComparison( sName ) )
	{
		cout << "ERROR: A comparison MUST have a unique name! this name is used more than once, " << sName << endl;
		return false;
	}
	
	// check if this indirect comparison uses two valid direct comparison; 
	ProteinDirectComparison tempDirectComparison0;
	ProteinDirectComparison tempDirectComparison1;
	if( !getDirectComparison( tempProteinIndirectComparison.getDirectComparisonName0(), tempDirectComparison0 ) )
	{
		cout << "ERROR: indirect comparison has an invalid DIRECT_COMPARISON_NAME_A " 
			<< tempProteinIndirectComparison.getDirectComparisonName0() << endl;
		return false;
	}
	
	if( !getDirectComparison( tempProteinIndirectComparison.getDirectComparisonName1(),  tempDirectComparison1) )
	{
		cout << "ERROR: indirect comparison has an invalid DIRECT_COMPARISON_NAME_B " 
			<< tempProteinIndirectComparison.getDirectComparisonName1() << endl;
		return false;
	}
	if( !tempProteinIndirectComparison.testDirectComparison(tempDirectComparison0, tempDirectComparison1) )
	{
		return false;
	}
			
	vIndirectComparison.push_back( tempProteinIndirectComparison );
	return true;
}

const vector< ProteinDirectComparison > & ProteinCombined::getProteinDirectComparisonVector()
{
	return vDirectComparison;
}

const vector< ProteinIndirectComparison > & ProteinCombined::getProteinIndirectComparisonVector()
{
	return vIndirectComparison;
}

bool ProteinCombined::getDirectComparison( string sComparisonName, ProteinDirectComparison & comparison )
{
	int i = 0;
	for( i = 0; i < vDirectComparison.size(); ++i )
	{
		if( vDirectComparison[i].getName() == sComparisonName )
		{
			comparison = vDirectComparison[i];
			return true;
		}
	}
	return false;

}

bool ProteinCombined::getIndirectComparison( string sComparisonName, ProteinIndirectComparison & comparison )
{
	int i = 0;
	for( i = 0; i < vIndirectComparison.size(); ++i )
	{
		if( vIndirectComparison[i].getName() == sComparisonName )
		{
			comparison = vIndirectComparison[i];
			return true;
		}
	}
	return false;

}

bool ProteinCombined::IsThereDirectComparison( string sComparisonName )
{
	ProteinDirectComparison temp;
	return getDirectComparison( sComparisonName, temp );
}

bool ProteinCombined::IsThereIndirectComparison( string sComparisonName )
{
	ProteinIndirectComparison temp;
	return getIndirectComparison( sComparisonName, temp );
}

void ProteinCombined::computeValidity()
{
	// for a ProteinCombined to be valid, it has to have at least a valid direct comparison or indirect comparison
	bValidity = false;
	int i = 0;
	for( i = 0; i < vDirectComparison.size(); ++i )
	{
		if( vDirectComparison[i].isValid() )
		{
			bValidity = true;
		}
	}
	for( i = 0; i < vIndirectComparison.size(); ++i )
	{
		if( vIndirectComparison[i].isValid() )
		{
			bValidity = true;
		}
	}
}

bool LessProteinCombined::operator() ( ProteinCombined * pProtein1, ProteinCombined * pProtein2 ) const
{
	if( sKey == "locus" )
	{
		if( pProtein1->getLocus() < pProtein2->getLocus() )
			return true;
		else
			return false;
	}
	else if ( sKey == "description" )
	{
		if( pProtein1->getDescription() < pProtein2->getDescription() )
			return true;
		else
			return false;
	}
	else
	{
		if( pProtein1->getLocus() < pProtein2->getLocus() )
			return true;
		else
			return false;
	}
}

