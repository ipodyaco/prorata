
#ifndef PRORATATANDEMMASSSPECTRUM_H
#define PRORATATANDEMMASSSPECTRUM_H

#include "proRataGraph.h"
#include "qwt_data.h"
#include "qwt_plot_picker.h"
#include "qwt_plot_zoomer.h"
#include "qwt_legend.h"
#include "qwt_legend_item.h"
#include "peptideRatio.h"
#include "directoryStructure.h"
#include "proRataMassSpectrumData.h"
#include "mzReader.h"
#include <QString>
#include <QLabel>
#include "proRataMassSpectrum.h"
#include "isotopologue.h"

class ProRataTandemMassSpectrum : public ProRataMassSpectrum
{
	Q_OBJECT;
	public:
		ProRataTandemMassSpectrum( QWidget * qwParent = 0 );
		~ProRataTandemMassSpectrum();

	public slots:
		virtual void peptideUpdated( PeptideRatio * );
		virtual void updateMSGraph( long );
		virtual void cleanGraph();
	protected:
		string sSequence;
		int iChargeState;
		vector< unsigned long int > viMS2ScanNumber;
		vector< Isotopologue * > vpIsotopologue;

		void markIonSeries( vector< double > vdIonMZ, const QColor & qcolor, const QString & qsName );

};
#endif //PRORATATANDEMMASSSPECTRUM_H

