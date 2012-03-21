
#ifndef PRORATAMAINWIDGET_H
#define PRORATAMAINWIDGET_H

#include <QWidget>
#include <QDockWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QString>
#include <QSplitter>
#include <QFile>
#include <QDir>
#include <QIcon>

#include "proRataTablePane.h"
#include "proRataImagePane.h"
#include "proRataGraphPane.h"
#include "proRataTextPane.h"
#include "proRataTextArea.h"
#include "proRataXmlProcessor.h"


#include "proRataLikelihoodGraph.h"
#include "proRataChromatogramGraph.h"
#include "proRataMassSpectrum.h"
#include "proRataTandemMassSpectrum.h"
#include "proRataPCAGraph.h"
#include "proRataPPCGraph.h"
#include "proRataSequenceCoverageGraph.h"
#include "proRataConfig.h"
#include "proRataMassSpectrumData.h"

class ProRataMainWidget : public QWidget
{
	Q_OBJECT;
	public:
		ProRataMainWidget( const QString qsFname, QWidget * qwParent = 0 );
		ProRataMainWidget( QWidget * qwParent = 0 );
		~ProRataMainWidget();

		void setFilename( const QString &qsFname );
		void showTables( bool );
		void showText( bool );
		void showGraphs( bool );

		void updateUI();

		bool redirectSearchSlot( bool );

	public slots:
		void newSearchString( QString );

	signals:
		void updateStatus( int );
		void updateStatusMessage( QString );

	private:
		QWidget *qwTableTop;
		ProRataTablePane * prtpTables;
		ProRataGraphPane * prgpGraphs;
		ProRataTextPane * prtpText;
		ProRataSearchPane * prspFind;

		ProRataProteinTable *prtProtein;
		ProRataPeptideTable *prtPeptide;

		ProRataXmlProcessors prxp;
		ProteomeInfo *mainProteomeInfo;
		ProRataMassSpectrumData *mainMassSpecData;

		QString qsFilename;

};
#endif //PRORATAMAINWIDGET_H
