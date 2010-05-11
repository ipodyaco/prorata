
#ifndef PRORATACHROMATOGRAM_H
#define PRORATACHROMATOGRAM_H

#include "proRataGraph.h"
#include "peptideRatio.h"
#include "qwt_data.h"
#include "qwt_plot_picker.h"
#include "qwt_plot_marker.h"
#include "qwt_double_rect.h"
#include <qwt_legend.h>
#include <QString>
#include <QLabel>

#include <QWidget>
#include <QWheelEvent>

class ProRataChromatogramGraph : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataChromatogramGraph( QWidget * qwParent = 0  );
		~ProRataChromatogramGraph();

		void setStartBoundaryMarker( double dValue );
		void setEndBoundaryMarker( double dValue );
		void setMSMSScanMarker( double dValue );

		void setReferenceData( const QwtArray<double> &, 
				const QwtArray<double> & );
		void setTreatmentData( const QwtArray<double> &, 
				const QwtArray<double> & );
		
	signals:
		void MS1selected( long );
		void MS2selected( long );

	public slots:
		virtual void peptideUpdated( PeptideRatio * );
		virtual void proteinUpdated( ProteinRatio * );
		virtual void leftPressed( const QwtDoublePoint &pos );
		virtual void cleanGraph();

	protected:
		void wheelEvent(QWheelEvent *event);
	private:
		QwtArrayData *qwtdReferenceData;
		QwtArrayData *qwtdTreatmentData;
		
		QwtLegend * qwtlLegend;

		QwtPlotCurve *qwtpcReferenceCurve;
		QwtPlotCurve *qwtpcTreatmentCurve;

		QwtPlotMarker *qwtpmStartBoundary;
		QwtPlotMarker *qwtpmEndBoundary;

		unsigned long int iMS1ScanNumber;

		// variables for drawing the MS markers
		QwtPlotMarker *qwtpmFullScan;
		bool bFullScanMarkerExists;
		PeptideRatio *peppRatioCurrent;
		bool bChroExists;

};
#endif //PRORATACHROMATOGRAM_H
