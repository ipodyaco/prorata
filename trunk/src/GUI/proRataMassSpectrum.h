
#ifndef PRORATAMASSSPECTRUM_H
#define PRORATAMASSSPECTRUM_H

#include "proRataGraph.h"
#include "qwt_data.h"
#include "qwt_plot_picker.h"
#include "qwt_plot_zoomer.h"
#include "peptideRatio.h"
#include "directoryStructure.h"
#include "proRataMassSpectrumData.h"
#include "mzReader.h"
#include <QString>
#include <QLabel>
#include <qwt_legend.h>
//class ProRataMassSpectrumData;

class ProRataMassSpectrum : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataMassSpectrum( QWidget * qwParent = 0 );
		~ProRataMassSpectrum();
		void setData( ProRataMassSpectrumData * pMassSpecDataInput  );

	public slots:
		virtual void peptideUpdated( PeptideRatio * );
		virtual void updateMSGraph( long );
		virtual void cleanGraph();
	protected:
		
		void setFullScanData( const QwtArray<double> &, const QwtArray<double> & , const QString & qsName );
		void setMZRange(  vector< float > vfLower, vector< float > vfUpper, const QColor & qcolor, const QString & qsMZRangeName  );

		QwtArrayData *qwtdFullScanData;
		QwtPlotCurve *qwtpcFullScanCurve;

	//	QwtLegend * qwtlLegend;
		PeptideRatio *ppepActiveRatio;
		
		ProRataMassSpectrumData * pMassSpecData;

		bool bZoomBaseSet;

		QwtPlotZoomer *qwtZoomer;

};
#endif //PRORATAMASSSPECTRUM_H
