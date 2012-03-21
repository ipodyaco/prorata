
#ifndef PRORATAISOTOP_H
#define PRORATAISOTOP_H

#include <QVariant>
#include <QPushButton>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QTableWidget>
#include <QHeaderView>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QStringList>
#include <vector>
#include <QComboBox>

#include <QApplication>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QStandardItemModel>
#include <QTableView>
#include "delegate.h"

using namespace std;

class ProRataIsotopologue : public QWidget
{
	Q_OBJECT

	public:
		
		ProRataIsotopologue( QWidget* qwParent = 0, Qt::WFlags qwfFl = 0 );
		~ProRataIsotopologue();

		void buildUI();
		void setName( const QString & );
		void insertRow( QString qsRowName, vector< int > viRowContent );
		void getData( QStringList & qsRowNameListInput, vector< vector< int > > & vviTableContentInput );
		QString getName() { return qleName->text(); }
		
	//	void setData( int iColumn, const QStringList & );
	//	const QStringList & getValues();

	public slots:
		void insertPTMrow();
		
	protected:

		QString rowData( int iRow );
		
		QStringList qslRowNameList;

		QGroupBox* qgbIsotop;

		QLabel* qtlName;
		QLineEdit* qleName;
	//	QTableWidget* qtTable;
	
		QSpacerItem * qSpacer;
		// widgets for inserting rows
		QLabel* qtlInsertLabel;
		QComboBox* qlePTMletter;
		QPushButton* qpbAddPTM;
		

		// derived item model
		QStandardItemModel *model;
		QTableView *tableView;
		SpinBoxDelegate *delegate;

		QGridLayout* qglProRataIsotopologueLayout;
		QGridLayout* qgbIsotopLayout;

		int iColCt;
		int iRowCt;

		QStringList qslValues;


};
#endif // PRORATAISOTOPOLOGUE_H


