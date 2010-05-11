
#ifndef PRORATAPROTEINTABLE_H
#define PRORATAPROTEINTABLE_H

#include <QTableView>
#include <QVBoxLayout>
#include <QAbstractItemModel>
#include <QStandardItemModel>
#include <QModelIndex>

#include <math.h>

#include <vector>
#include <string>

#include "proRataTable.h"
#include "proRataXmlProcessor.h"
#include "proRataProteinTuple.h"

#include "proteinInfo.h"
using namespace std;

class ProRataProteinTable : public ProRataTable
{
	Q_OBJECT;
	
	public:
		ProRataProteinTable( ProteomeInfo *, QWidget * qwParent = 0  );
		~ProRataProteinTable();

		virtual void setXmlProcessor( ProRataXmlProcessors * );
		virtual void populateTable( const vector< ProteinInfo *> & );

	private slots:
		void rowClicked( const QModelIndex & );
		void colClickedForSorting( int iColumn );


	signals:
		void proteinClicked( const QString & );
		void proteinClicked( ProteinInfo * );
		void proteinClicked( ProteinRatio * );
		void flushGraph();


	private:

		int iCurrentColumnSorted;
		bool bIsCurrentSortAscending;
		void populateTable();
		void sortColumn( int iColumn, bool bIsAscending );

		ProRataXmlProcessors *prxpProcessor;
		ProRataProteinTuple *prProteinTuple;
		int iRowCount;
		ProteinRatio * pproRatio;

		vector< ProteinInfo *> vproInfo;


};
#endif //PRORATAPROTEINTABLE_H

