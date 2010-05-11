
#ifndef PRORATATABLE_H
#define PRORATATABLE_H

#include <QTableView>
#include <QVBoxLayout>
#include <QAbstractItemModel>
#include <QStandardItemModel>
#include <QModelIndex>
#include <QStringList>
#include <QModelIndex>
#include <QLabel>
#include <QHeaderView>

#include "titleLabel.h"
#include "proRataXmlProcessor.h"
#include "proteomeInfo.h"

#include <math.h>

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class ProRataTable : public QWidget
{
	Q_OBJECT;
	
	public:

		ProRataTable( ProteomeInfo *, QWidget * qwParent = 0  );
		~ProRataTable();

		void setupModel( const QStringList & );
		void setTableTitle( const QString & qsTitle )
		{	qlTableTitle->setText( qsTitle ); }

		virtual void setXmlProcessor( ProRataXmlProcessors * ){}
		virtual void setIndicator( int );

		virtual void toggleDescending()
		{	if (bDescending)
					bDescending = false;
			else
					bDescending = true;
		}

		virtual bool isDescending()
		{	return bDescending;	}

	public slots:
		void colClickedForSorting( int iColumn );

	protected:

		ProteomeInfo *mainProteomeInfoInstance;
		QTableView *qtwTable;
		QAbstractItemModel *qaiModel;
		
	private:

		TitleLabel *qlTableTitle;
		QString qsTableTitle;
		bool bDescending;

};
#endif //PRORATATABLE_H

