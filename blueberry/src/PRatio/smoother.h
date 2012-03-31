
#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

extern "C"
{
	extern void savgol( float *, int, int, int, int, int );
}

class Smoother
{
	public:
		Smoother( int iWinSize, int iOrder );
		~Smoother();
		void setWindowSize( int iWinSize ) { 
			iWindowSize = iWinSize;}

		void setDesiredOrder( int iOrder ) {
			iDesiredOrder = iOrder;}

		int smoothen( std::vector<double> & );

	private:

		int iWindowSize;
		int iDesiredOrder;


};

#endif //SMOOTHER_H
