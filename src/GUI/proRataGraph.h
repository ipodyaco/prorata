
#ifndef PRORATAGRAPH_H
#define PRORATAGRAPH_H

#include "qwt_plot.h"
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_picker.h>
#include <qwt_symbol.h>
#include <qwt_legend.h>
#include <qwt_data.h>

#include <QPixmap>

#include "proRataExportImage.h"


#include "proteinRatio.h"
#include "peptideRatio.h"

class ProRataGraph : public QwtPlot
{
	Q_OBJECT;
	public:
		ProRataGraph( QWidget * qwParent = 0  );
		virtual ~ProRataGraph();
		void saveToFile(const QString &qsFilename, const char *cpFormat);
		void setPosition( int iPos = 0 )
		{	iPosition = iPos;	}
		int getPosition()
		{	return iPosition;	}

		QwtPlotPicker * getPicker()
		{	return qwtppPicker;	}

/*
		enum Type
		{
			LIKELIHOOD
		}
*/

	signals:
		void updated(QString, int);

	public slots:
		virtual void proteinUpdated( ProteinRatio * );
		virtual void peptideUpdated( PeptideRatio * );

		virtual void cleanGraph();

	protected:
	//	void legendStatus( bool );
	//	QwtLegend * qwtlLegend;
		QwtPlotPicker * qwtppPicker;
		int iPosition;

};
#endif //PRORATAGRAPH_H
