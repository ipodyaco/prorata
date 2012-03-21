
#include "proRataProteinTable.h"
#include <QMessageBox>

ProRataProteinTable::ProRataProteinTable( ProteomeInfo * ptmInfoInst, QWidget * qwParent ) 
	: ProRataTable( ptmInfoInst, qwParent )
{
	iRowCount = 0;
	iCurrentColumnSorted = 0;
	bIsCurrentSortAscending = true;
	qtwTable->horizontalHeader()->setSortIndicatorShown( true );
	//qtwTable->horizontalHeader()->setSortIndicator( 0, Qt::AscendingOrder );
	qtwTable->horizontalHeader()->setSortIndicator( iCurrentColumnSorted, Qt::AscendingOrder );
	qtwTable->horizontalHeader()->resizeSection( 0, 80 );
	qtwTable->horizontalHeader()->resizeSection( 1, 60 );
	qtwTable->horizontalHeader()->resizeSection( 2, 60 );

	connect( qtwTable, SIGNAL( clicked( const QModelIndex & ) ), 
			this, SLOT( rowClicked( const QModelIndex & ) ) );

	connect( qtwTable->horizontalHeader(), SIGNAL( sectionClicked( int ) ),
		this, SLOT( colClickedForSorting( int ) ) );



	pproRatio = NULL;
}

ProRataProteinTable::~ProRataProteinTable()
{
}

void ProRataProteinTable::setXmlProcessor( ProRataXmlProcessors * prxpProc )
{


}

void ProRataProteinTable::rowClicked( const QModelIndex & qmiClickedItem )
{
	int iRow = qmiClickedItem.row();

	if ( iRow < 0 )
	{
		return;
	}

	qtwTable->selectRow( iRow );

	delete pproRatio;

	pproRatio = new ProteinRatio;

	vproInfo.at(iRow)->setProteinRatio( pproRatio );

	//emit proteinClicked( qsTempLocus );
	emit proteinClicked( vproInfo.at( iRow ) );
	emit proteinClicked( pproRatio );

}

void ProRataProteinTable::populateTable( const vector< ProteinInfo *> &vpProteinInfo )
{
	//emit flushGraph();

	vproInfo.clear();
	vproInfo.resize( vpProteinInfo.size() );
	copy( vpProteinInfo.begin(), vpProteinInfo.end(), vproInfo.begin() );
	sortColumn( iCurrentColumnSorted, bIsCurrentSortAscending );
	//populateTable();
	//qtwTable->resizeColumnToContents( 0 );
	//qtwTable->resizeColumnToContents( 1 );
	//qtwTable->resizeColumnToContents( 2 );
	//qtwTable->resizeColumnToContents( 3 );
	//qtwTable->horizontalHeader()->resizeSection( 0, 70 );
	qtwTable->horizontalHeader()->setStretchLastSection( true );


}

void ProRataProteinTable::populateTable()
{
	qaiModel->removeRows( 0, qaiModel->rowCount( QModelIndex() ), QModelIndex() );
	iRowCount = 0;

	int iTupleSize = vproInfo.size();
	
//	qaiModel->beginInsertRows();
	qaiModel->insertRows( iRowCount, iTupleSize, QModelIndex() );
//	qaiModel->endInsertRows();
	
	for( int i = 0; i < iTupleSize; i++ )
	{
		// Add Locus
		qaiModel->setData( qaiModel->index( iRowCount, 0, QModelIndex() ),
				vproInfo.at(i)->getLocus().c_str() );

		// Add Log2Ratio
		qaiModel->setData( qaiModel->index( iRowCount, 1, QModelIndex() ),
				QString::number( vproInfo.at(i)->getLog2Ratio() ) );

		// Add Confidence Interval Width
		qaiModel->setData( qaiModel->index( iRowCount, 2, QModelIndex() ),
				QString::number( vproInfo.at(i)->getUpperLimitCI() - 
				vproInfo.at(i)->getLowerLimitCI() ) );

		// Add Description
		qaiModel->setData( qaiModel->index( iRowCount, 3, QModelIndex() ),
				vproInfo.at(i)->getDescription().c_str() );

		qtwTable->resizeRowToContents( iRowCount );
		iRowCount++;
	}
}

void ProRataProteinTable::colClickedForSorting( int iColumn )
{
	if( iColumn == iCurrentColumnSorted )
	{
		if( bIsCurrentSortAscending )
		{
			sortColumn( iColumn, false );

		}
		else
		{
			sortColumn( iColumn, true );
		}
	}
	else
	{
		sortColumn( iColumn, true );
	}
	
}


void ProRataProteinTable::sortColumn( int iColumn, bool bIsAscending )
{
	if( bIsAscending )
	{
		iCurrentColumnSorted = iColumn;
		bIsCurrentSortAscending = true;
		qtwTable->horizontalHeader()->setSortIndicator( iColumn, Qt::AscendingOrder );
		switch ( iColumn )
		{
		case (0):
			mainProteomeInfoInstance->sortProteinInfo( vproInfo, "locus" );
			break;
		case (1):
			mainProteomeInfoInstance->sortProteinInfo( vproInfo, "log2Ratio" );
			break;
		case (2):
			mainProteomeInfoInstance->sortProteinInfo( vproInfo, "widthCI" );
			break;
		case (3):
			mainProteomeInfoInstance->sortProteinInfo( vproInfo, "description" );
			break;
		default:
			break;
		}

	}
	else
	{
		iCurrentColumnSorted = iColumn;
		bIsCurrentSortAscending = false;
		qtwTable->horizontalHeader()->setSortIndicator( iColumn, Qt::DescendingOrder );
		switch ( iColumn )
		{
		case (0):
			mainProteomeInfoInstance->sortProteinInfoDescending( vproInfo, "locus" );
			break;
		case (1):
			mainProteomeInfoInstance->sortProteinInfoDescending( vproInfo, "log2Ratio" );
			break;
		case (2):
			mainProteomeInfoInstance->sortProteinInfoDescending( vproInfo, "widthCI" );
			break;
		case (3):
			mainProteomeInfoInstance->sortProteinInfoDescending( vproInfo, "description" );
			break;
		default:
			break;
		}

	}
	emit flushGraph();
	populateTable();


}

