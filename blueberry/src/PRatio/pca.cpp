
#include "pca.h"

PCA::PCA( int iRow, int iColumn )
{
	setRow( iRow );
	setColumn( iColumn );
}

PCA::~PCA()
{

}

int PCA::calculate( const std::vector<double> &vdColumnOne, const std::vector<double> &vdColumnTwo )
{
	int i, j;
	float **data, **symmat, *evals, *interm;

	data = pcaMatrix( iRow, iColumn );


	for (i = 1; i <= iRow; i++)
	{
		data[i][1] = vdColumnOne.at(i-1);
		data[i][2] = vdColumnTwo.at(i-1);
	}

	symmat = pcaMatrix(iColumn, iColumn);

	covcol(data, iRow, iColumn, symmat);


	evals = pcaVector(iColumn); 
	interm = pcaVector(iColumn);

	tred2(symmat, iColumn, evals, interm);

	tqli(evals, interm, iColumn, symmat);


	iEigenVectorCount = iEigenValueCount = iColumn;


	for (j = iColumn; j >= 1; j--) 
	{
		vdEigenValues.push_back( evals[j] );
	}

	for (j = 1; j <= iColumn; j++) 
	{
		for (i = 1; i <= iColumn; i++)  
		{
			vdEigenVectors.push_back( symmat[j][iColumn - i + 1] );
		}
	}

	free_pcaMatrix(data, iRow, iColumn);
	free_pcaMatrix(symmat, iColumn, iColumn);
	free_pcaVector(evals, iColumn);
	free_pcaVector(interm, iColumn);

	return 0;
}
const vector<double> & PCA::getEigenValues()
{
	return vdEigenValues;
}

const vector<double> &  PCA::getEigenVectors()
{
	return vdEigenVectors;
}

