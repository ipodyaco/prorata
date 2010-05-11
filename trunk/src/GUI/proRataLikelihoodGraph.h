
#ifndef PRORATALIKELIHOODGRAPH_H
#define PRORATALIKELIHOODGRAPH_H

#include "proRataGraph.h"
#include "qwt_data.h"
#include <QString>
#include <QLabel>

class ProRataLikelihoodGraph : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataLikelihoodGraph( QWidget * qwParent = 0  );
		~ProRataLikelihoodGraph();

		void setLowerCI( double dValue );
		void setUpperCI( double dValue );
		void setCutoff( double dCutoff );
		/*void setType( Type gType = LIKELIHOOD)
		{	graphType = Type;	}
		*/

		/*
		void setLikelihoodData( const double *, const double *, int );
		void setPeptidePoints( const double *, const double *, int );
		*/

		void setLikelihoodData( const QwtArray<double> &, 
				const QwtArray<double> & );
		void setPeptidePoints( const QwtArray<double> &, 
				const QwtArray<double> & );


	signals:
	//	void rectMLESelected( const QwtDoubleRect );
		void pointMLESelected( const QwtDoublePoint );
		//void updated();
		
	public slots:
		virtual void proteinUpdated( ProteinRatio * );
		virtual void peptideUpdated( PeptideInfo * pPepInfo );

	private slots:	
	//	virtual void relayRectSelected( const QwtDoubleRect & );
		virtual void relayPointSelected( const QwtDoublePoint & );
		virtual void cleanGraph();

	
		
	private:
		QwtArrayData *qwtdLikelihoodData;
		QwtArrayData *qwtdPeptidePointsData;

		QwtPlotCurve *qwtpcLikelihoodCurve;
		QwtPlotCurve *qwtpcPeptidePointsCurve;

		QwtPlotCurve *qwtpcSelectedPeptide;

		QwtPlotMarker *qwtpmLowerCIMarker;
		QwtPlotMarker *qwtpmUpperCIMarker;
		QwtPlotMarker *qwtpmLikelihoodMarker;

		QwtSymbol qwtsPeptidePointsSym;
		QwtSymbol qwtsPeptideSelectedSym;

		ProteinRatio *propRatioCurrent;
		bool bPreviousSelected;
		bool bPopulated;


};
#endif //PRORATALIKELIHOODGRAPH_H
