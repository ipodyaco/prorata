
#include "smoother.h"

#include <iostream>

Smoother::Smoother( int iWinSize, int iOrder )
{
	setWindowSize( iWinSize );
	setDesiredOrder( iOrder );
}

Smoother::~Smoother()
{

}

int Smoother::smoothen( std::vector<double> & vdDataColumn )
{
	int iIndex;
	double dSum = 0.0;
	int iDataSize = vdDataColumn.size();
	int i, j, iInputIndex, iOutputIndex, iCoeffIndex;
	std::vector<double> vdOutput( iDataSize );

	// Need some error catching mechanism
	float *fpCoeffArray = new float[ iWindowSize + 1 ];
	float *fpCoeffArrayOrdered = new float[ iWindowSize ];

	int iDerivativeOrder = 0;
	int iLeftDataPoints = ( iWindowSize - 1 ) / 2 ;
	int iRightDataPoints = ( iWindowSize - 1 ) / 2 ;
	int iTotalPoints = iWindowSize;

	savgol( fpCoeffArray, iTotalPoints, iLeftDataPoints, 
			iRightDataPoints, iDerivativeOrder, iDesiredOrder );


	for (  i = 0; i < iWindowSize ; i++ )
	{
		iIndex = ( iTotalPoints / 2 ) + 1 - i ;
		if ( iIndex <= 0 )
			iIndex += iTotalPoints;
		fpCoeffArrayOrdered[i] = fpCoeffArray[iIndex];
	}

	

	for ( i = 0; i < iDataSize; i++ )
	{
		dSum = 0.0;

		for ( j = 0; j < iWindowSize; j++ )
		{
			iInputIndex = ( i + j ) % iDataSize ;
			iCoeffIndex = j;
			dSum = dSum + fpCoeffArrayOrdered[iCoeffIndex] * 
				vdDataColumn[iInputIndex];

		}
		iOutputIndex = ( i + (  iWindowSize / 2 )  ) % iDataSize;
		vdOutput[iOutputIndex] = dSum;
	}

	std::copy( vdOutput.begin(), vdOutput.end(), vdDataColumn.begin() );


	// Need some error catching mechanism
	delete [] fpCoeffArray;
	delete [] fpCoeffArrayOrdered;

	return 0;
}
