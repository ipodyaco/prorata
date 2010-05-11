
#ifndef PRORATATABLEPANE_H
#define PRORATATABLEPANE_H

#include "proRataProteinTable.h"
#include "proRataPeptideTable.h"
#include "proRataXmlProcessor.h"
#include "proRataSearchPane.h"

#include <QSplitter>

class ProRataTablePane : public QWidget
{
	Q_OBJECT;

	public:
		ProRataTablePane( ProRataXmlProcessors * prxpProcessor, QWidget * qwParent = 0  );
		ProRataTablePane( QWidget * qwParent = 0  );
		~ProRataTablePane();

		void setProteinTable( ProRataTable * prtTable );
		void setPeptideTable( ProRataTable * prtTable );

/*
	public slots:
		void hideSearchPane(bool);
		void reEmitSlot( QString );


	signals:
		void reEmitSearchString( QString  );
*/

	private:

		void buildUI();
    	QVBoxLayout *qvbLayout;
		ProRataProteinTable *prtProtein;
		ProRataPeptideTable *prtPeptide;
		ProRataSearchPane *prspFind;
		ProteomeInfo *mainProteomeInfoInstance;
		QSplitter *qsTableDivider;

};

#endif //PRORATATABLEPANE_H
