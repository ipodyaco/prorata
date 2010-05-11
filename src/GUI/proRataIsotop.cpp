#include "proRataIsotop.h"


ProRataIsotopologue::ProRataIsotopologue( QWidget* qwParent, Qt::WFlags qwfFl )
	: QWidget( qwParent, qwfFl )
{
	iColCt = 12;
	iRowCt = 0;

	buildUI();
}

ProRataIsotopologue::~ProRataIsotopologue()
{
}

void ProRataIsotopologue::buildUI()
{
	qglProRataIsotopologueLayout = new QGridLayout;

	qgbIsotop = new QGroupBox;
	qgbIsotopLayout = new QGridLayout;
	qgbIsotopLayout->setAlignment( Qt::AlignTop );
	qgbIsotop->setLayout( qgbIsotopLayout );

	qtlName = new QLabel( qgbIsotop );
	qgbIsotopLayout->addWidget( qtlName, 0, 0 );
	qtlName->setText( tr( "Name:" ) );
	qleName = new QLineEdit( qgbIsotop );
	qgbIsotopLayout->addWidget( qleName, 0, 1 );

	qSpacer = new QSpacerItem(25, 0, QSizePolicy::Expanding,QSizePolicy::Minimum);
	qgbIsotopLayout->addItem( qSpacer, 0, 2 );
	
	qtlInsertLabel = new QLabel( qgbIsotop );
	qgbIsotopLayout->addWidget( qtlInsertLabel, 0, 3 );
	qtlInsertLabel->setText( tr( " Add PTM:" ));

	qlePTMletter = new QComboBox( qgbIsotop );
	qgbIsotopLayout->addWidget( qlePTMletter, 0, 4 );
	qlePTMletter->insertItem( 0, "   ");
	qlePTMletter->insertItem( 1, "*");
	qlePTMletter->insertItem( 2, "@");
	qlePTMletter->insertItem( 3, "#");
	qlePTMletter->insertItem( 4, "&");
	qlePTMletter->insertItem( 5, "?");
	qlePTMletter->insertItem( 6, ">");
	qlePTMletter->insertItem( 7, "<");
	qlePTMletter->insertItem( 8, "%");
	qlePTMletter->insertItem( 9, "!");
	qlePTMletter->insertItem( 10, "$");
	qlePTMletter->insertItem( 11, "=");
	qlePTMletter->insertItem( 12, "(");
	qlePTMletter->insertItem( 13, ")");
	qlePTMletter->insertItem( 14, "{");
	qlePTMletter->insertItem( 15, "}");
	qlePTMletter->insertItem( 16, "[");
	qlePTMletter->insertItem( 17, "]");
	qlePTMletter->insertItem( 18, "+");
	qlePTMletter->insertItem( 19, "\"");
	qlePTMletter->insertItem( 20, "~");
	qlePTMletter->insertItem( 21, ",");
	qlePTMletter->insertItem( 22, "/");
	qlePTMletter->insertItem( 23, "|");
//	qlePTMletter->setMaxLength( 1 );

	qpbAddPTM = new QPushButton( qgbIsotop );
	qpbAddPTM->setText( tr( "Insert Row" ));
	qgbIsotopLayout->addWidget( qpbAddPTM, 0, 5 );

	connect( qpbAddPTM, SIGNAL( clicked() ), this, SLOT( insertPTMrow() ) );

	/*
	qtTable = new QTableWidget( qgbIsotop );
	qtTable->setColumnCount( iColCt );
	qtTable->setRowCount( iRowCt );
	qtTable->setEditTriggers( QAbstractItemView::AllEditTriggers );
	qtTable->setTabKeyNavigation( true );
	qtTable->setAlternatingRowColors( true);
//	qtTable->horizontalHeader()->setStretchLastSection( true );

	int i;
	int j;

	// make all cell aligned in the center
	for( i = 0; i < iRowCt; ++i )
	{
		for( j = 0; j < iColCt; ++j )
		{
		//	qtTable->item( i, j )->setTextAlignment( Qt::AlignCenter );
		}
	}

	// to make the first column of residue name not editable
	for( i = 0; i < iRowCt - 3; ++i )
	{
	//	qtTable->item( i, 0 )->setFlags( Qt::ItemIsSelectable ); 
	}
	*/

	model = new QStandardItemModel(0, iColCt);

	int i;
	QStringList qslHeaders;
	qslHeaders << 
		"C12" << "H1 " << "O16" << "N14" << "P31" << "S32" <<
		"C13" << "H2 " << "O18" << "N15" << "P32" << "S34";
	for (i = 0; i < qslHeaders.size(); ++i)
	{
		model->setHeaderData( i, Qt::Horizontal, 
				tr(qslHeaders.at( i ).toAscii() ) );
	}

	tableView = new QTableView;
	tableView->setModel(model);

	
	delegate = new SpinBoxDelegate;
	tableView->setItemDelegate(delegate);
	tableView->setEditTriggers( QAbstractItemView::AllEditTriggers );
	tableView->setTabKeyNavigation( true );
	tableView->setAlternatingRowColors( true);
	
	/*
	for (int row = 0; row < 25; ++row) {
		for (int column = 0; column < 12; ++column) {
		    QModelIndex index = model->index(row, column, QModelIndex());
		    model->setData(index, QVariant((row+1) * (column+1)));
		}
	}
	*/
	
	qgbIsotopLayout->addWidget( tableView, 1, 0, 1, 6 );
//	qgbIsotopLayout->addWidget( qtTable, 1, 0, 1, 2 );

	qglProRataIsotopologueLayout->addWidget( qgbIsotop, 0, 0 );

	qgbIsotop->setTitle( tr( "Residue Atomic Composition Table" ) );
	for( i = 0; i < 12; i++ )
		tableView->resizeColumnToContents( i );
	
	setLayout( qglProRataIsotopologueLayout );
}


void ProRataIsotopologue::setName( const QString & qsName )
{
	qleName->setText( qsName );
}

void ProRataIsotopologue::insertRow( QString qsRowName, vector< int > viRowContent )
{
	qslRowNameList.append( qsRowName );
	
	model->insertRows( model->rowCount(), 1 );
	model->setHeaderData(  ( model->rowCount() - 1 ), Qt::Vertical, tr( qsRowName.toAscii() ) );
	for (int column = 0; column < iColCt; ++column) {
	    QModelIndex index = model->index( ( model->rowCount() - 1 ), column, QModelIndex());
	    model->setData(index, QVariant( viRowContent[column] ));
	}	
	tableView->resizeRowToContents( model->rowCount() - 1 );	

}

void ProRataIsotopologue::insertPTMrow()
{
	QString qsName = qlePTMletter->currentText();
	if( qsName == "" || qsName == "   " )
		return;

	qlePTMletter->removeItem( qlePTMletter->currentIndex() );
	qlePTMletter->setCurrentIndex( 0 );
	
	vector< int > viInitialRow;
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	viInitialRow.push_back( 0 );
	insertRow( qsName, viInitialRow );

}

void ProRataIsotopologue::getData( QStringList & qsRowNameListInput, vector< vector< int > > & vviTableContentInput )
{
	delegate->commitCurrentEditor( model );

//	tableView->edit( tableView->currentIndex() );
//	tableView->closePersistentEditor( tableView->currentIndex() );
//	QMessageBox::information(this, "after ", QString::number( tableView->state() ) );
	qsRowNameListInput.clear();
	vviTableContentInput.clear();
	qsRowNameListInput = qslRowNameList;
	
	int i;
	int j;
	for( i = 0; i < qslRowNameList.size(); i++ )
	{
		vector< int > viRowContent;
		for( j = 0; j < iColCt ; j++ )
		{
			QModelIndex index = model->index( i, j, QModelIndex());
			QVariant qvCell = model->data(index);
			viRowContent.push_back( qvCell.toInt() );
		}
		vviTableContentInput.push_back( viRowContent );

	}
	
}

/*
void ProRataIsotopologue::setData( int iColumn, const QStringList &qslData )
{

	for( int i = 0; i < iRowCt - 3; i++ )
	{
		QTableWidgetItem* tempItem = new QTableWidgetItem(  qslData.at(i) );
		tempItem->setTextAlignment( Qt::AlignCenter );
		qtTable->setItem( i, iColumn, tempItem );
		//qtTable->setItem( i, iColumn, new QTableWidgetItem( qslData.at(i) ) );
	}

	for(int i = iRowCt - 3; i < iRowCt; i++ )
	{
		QTableWidgetItem* tempItem = new QTableWidgetItem(  QString("") );
		tempItem->setTextAlignment( Qt::AlignCenter );
		qtTable->setItem( i, iColumn, tempItem );
		//qtTable->setItem( i, iColumn, new QTableWidgetItem( qslData.at(i) ) );
	}
}

const QStringList &  ProRataIsotopologue::getValues()
{
	qslValues.clear();

	qslValues << qleName->text();

	for( int i = 0; (i < iRowCt) && ( qtTable->item( i, 0 ) ); i++ )
	{
		qslValues << rowData( i );
	}

	return qslValues;
}

QString ProRataIsotopologue::rowData( int iRow )
{
	QString qsRowData;

	qsRowData = qtTable->item( iRow, 0 )->text();

	for( int i = 1; i < iColCt; i++ )
	{
		QString qsCellData = qtTable->item( iRow, i )->text();
		if (qsCellData == "" )
		{
			qsRowData = qsRowData + ",\t" + QString( "0" );
		}
		else
		{
			qsRowData = qsRowData + ",\t" + qsCellData;
		}
	}

	return qsRowData;
}
*/


