
#ifndef PRORATAPCAGRAPH_H
#define PRORATAPCAGRAPH_H

#include "proRataGraph.h"
#include "qwt_data.h"
#include "peptideRatio.h"
#include <QString>
#include <QLabel>

class ProRataPCAGraph : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataPCAGraph( QWidget * qwParent = 0  );
		~ProRataPCAGraph();

		void setPeptideData( const QwtArray<double> &, 
				const QwtArray<double> & );

		void setPrincipalComponent( double dX1, double dY1,
				double dX2, double dY2 );

		/*
		void setFirstPrincipalComponent( double dX1, double dY1,
				double slope );
		*/

	//	void setPrincipalComponents( double dAbundanceRatio, double dEigenValueRatio );

	public slots:
		virtual void peptideUpdated( PeptideRatio * );
		virtual void MS1selected( long );
		virtual void cleanGraph();
	private:

		bool rangeOfData( double * yMin, double * yMax, double *xMin, 
				double *xMax, double *xAvg, double *yAvg );

		QwtArrayData *qwtdPeptideData;
		QwtPlotCurve *qwtpcPeptideCurve;

		QwtSymbol qwtsPeptideSym;
		bool bSecondPCA;

		PeptideRatio *peppRatioCurrent;
		QwtPlotCurve *qwtpcSelectedMS1;
		QwtSymbol qwtsMS1Sym;
		bool bIsMS1Selected;

};
#endif //PRORATAPCAGRAPH_H
