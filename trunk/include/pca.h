
/*
 * Wrapper class internally uses Numerical receipes PCA algorithm
 */

#ifndef PCA_H
#define PCA_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

/*
 * C prototypes.
 * Should document later.
 */

extern "C"
{
	extern float **pcaMatrix( int, int );
	extern float *pcaVector( int );
	extern void free_pcaMatrix( float **, int, int );
	extern void free_pcaVector( float *, int );
	extern void covcol( float **, int, int, float ** );
	extern void tred2( float **, int, float *, float * );
	extern void tqli( float *, float *, int, float ** );
}


class PCA
{
	public:
		/*
		 * Create an instance by providing row and column count
		 */
		PCA( int iRow, int iColumn );
		~PCA();

		/*
		 * Utility functions setting Row/Column counts.
		 */
		void setRow( int iRowSize ) { iRow = iRowSize; }
		void setColumn( int iColumnSize ) { iColumn = iColumnSize; }

		/*
		 * Provide the matix and start calculate Eigen values
		 * and Eigen vectors.
		 */
		int calculate( const std::vector<double> &, const std::vector<double> & );


		int getNumberOfEigenValues() { return iEigenValueCount; }
		int getNumberOfEigenVectos() { return iEigenVectorCount; }

		/*
		 * Methods for geting Eigen values and Eigen vectors.
		 * Note the vectors have to be allocated before calling
		 * these methods. The size of the vectors can be got by
		 * above two methods.
		 */

		const vector<double> & getEigenValues();
		const vector<double> & getEigenVectors();


	private:
		int iRow;
		int iColumn;
		int iEigenValueCount;
		int iEigenVectorCount;
		vector<double> vdEigenValues;
		vector<double> vdEigenVectors;

};

#endif //PCA_H
