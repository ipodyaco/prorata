
#ifndef PRORATATEXTPANE_H
#define PRORATATEXTPANE_H

#include "proRataTextArea.h"
#include "proRataProteinTuple.h"
#include "proRataPeptideTuple.h"
#include "proRataProteinTable.h"
#include "proRataPeptideTable.h"

#include <QWidget>
#include <QString>
#include <QHBoxLayout>

class ProRataTextPane : public QWidget
{
	Q_OBJECT;

	public:
		ProRataTextPane( QWidget * qwParent = 0 );
		~ProRataTextPane();
		
		void addTextArea( ProRataTextArea * );
		void setProteinTable( ProRataProteinTable * );
		void setPeptideTable( ProRataPeptideTable * );

	private:

		ProRataTextArea *prtaText;
		QHBoxLayout *qhblMainLayout;
};
#endif //PRORATATEXTPANE_H

