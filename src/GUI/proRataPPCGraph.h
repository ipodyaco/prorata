
#ifndef PRORATAPPCGRAPH_H
#define PRORATAPPCGRAPH_H

#include "proRataGraph.h"
#include "qwt_data.h"
#include "peptideRatio.h"
#include <QString>
#include <QLabel>

class ProRataPPCGraph : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataPPCGraph( QWidget * qwParent = 0  );
		~ProRataPPCGraph();

		void setStartBoundaryMarker( double dValue );
		void setEndBoundaryMarker( double dValue );
		void setMSMSScanMarker( double dValue );

		void setPPCData( const QwtArray<double> &, 
				const QwtArray<double> & );

	public slots:
		virtual void peptideUpdated( PeptideRatio * );
		virtual void cleanGraph();
	private:
		QwtArrayData *qwtdPPCData;

		QwtPlotCurve *qwtpcPPCCurve;

		QwtPlotMarker *qwtpmStartBoundary;
		QwtPlotMarker *qwtpmEndBoundary;

};
#endif //PRORATAPPCGRAPH_H
