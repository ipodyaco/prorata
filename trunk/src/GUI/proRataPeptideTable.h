

#ifndef PRORATAPEPTIDETABLE_H
#define PRORATAPEPTIDETABLE_H

#include <QTableView>
#include <QVBoxLayout>
#include <QAbstractItemModel>
#include <QStandardItemModel>
#include <QModelIndex>

#include <math.h>
#include "qwt_double_rect.h"
#include <vector>
#include <string>

#include "proRataTable.h"
#include "proRataXmlProcessor.h"

using namespace std;

class ProRataPeptideTable : public ProRataTable
{
	Q_OBJECT;
	
	public:
		ProRataPeptideTable( ProteomeInfo *, QWidget * qwParent = 0  );
		~ProRataPeptideTable();

		virtual void setXmlProcessor( ProRataXmlProcessors * );
	//	virtual void populateTable( const vector< PeptideInfo *> & );
		
	public slots:
	//	void rectMLESelected( const QwtDoubleRect );
		void pointMLESelected( const QwtDoublePoint );

	private slots:
		void rowClicked( const QModelIndex & );
	//	void newProteinClicked( const QString & );
		void newProteinClicked( ProteinInfo * );
		void colClickedForSorting( int iColumn );
		void cleanUp();

	signals:
		void peptideClicked( PeptideRatio *);
		void peptideClicked( PeptideInfo *);
		void flushGraph();

	protected:
		void keyPressEvent ( QKeyEvent * event );


	private:
		void populateTable();
		void sortColumn( int iColumn, bool bIsAscending );
		ProRataXmlProcessors *prxpProcessor;
		ProRataPeptideTuple *prPeptideTuple;
		int iRowCount;
		PeptideRatio *ppepRatio;
		int iCurrentColumnSorted; 
		bool bIsCurrentSortAscending;
		vector<PeptideInfo *> vpepInfo;

};
#endif //PRORATAPEPTIDETABLE_H

