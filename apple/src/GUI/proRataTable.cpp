
#include "proRataTable.h"
#include <QMessageBox>

ProRataTable::ProRataTable( ProteomeInfo * ptmInfoInst, QWidget * qwParent ) : QWidget( qwParent )
{
	mainProteomeInfoInstance = ptmInfoInst;
	qlTableTitle = new TitleLabel( tr( "Proteomics Data" ) );
	qlTableTitle->tableHeader();

	qtwTable = new QTableView;
	qtwTable->horizontalHeader()->setStretchLastSection( true );
	
	qtwTable->setEditTriggers( QAbstractItemView::NoEditTriggers );
	qtwTable->setSelectionMode( QAbstractItemView::SingleSelection );
	qtwTable->setSelectionBehavior( QAbstractItemView::SelectRows );

	qaiModel = new QStandardItemModel( 10, 4, this );
	qtwTable->horizontalHeader()->setClickable( true );

	/*
	connect( qtwTable, SIGNAL( clicked( const QModelIndex & ) ), 
			this, SLOT( rowClicked( const QModelIndex & ) ) );
	*/

	qaiModel->removeRows( 0, qaiModel->rowCount( QModelIndex() ), QModelIndex() );

	qtwTable->setModel( qaiModel );
	//qtwTable->sortByColumn( 0 );

	QVBoxLayout *qvbLayout = new QVBoxLayout;
	qvbLayout->addWidget( qlTableTitle );
	qvbLayout->addWidget( qtwTable );

	qvbLayout->setSpacing( 0 );
	qvbLayout->setMargin( 0 );

	bDescending = false;

	//qvbLayout->addStretch( 1 );

	setLayout( qvbLayout );
}

ProRataTable::~ProRataTable()
{
}

void ProRataTable::setupModel( const QStringList & qslTableHeaders )
{
	for (int i = 0; i < qslTableHeaders.size(); ++i)
	{
		qaiModel->setHeaderData( i, Qt::Horizontal, 
				tr(qslTableHeaders.at( i ).toAscii() ) );
	}
}

void ProRataTable::colClickedForSorting( int iColumn )
{
}

void ProRataTable::setIndicator( int iColumn )
{
	if ( isDescending() )
	{
		qtwTable->horizontalHeader()->setSortIndicator( iColumn, Qt::DescendingOrder );
	}
	else
	{
		qtwTable->horizontalHeader()->setSortIndicator( iColumn, Qt::AscendingOrder );
	}
}
