
#include "proRataPeptideTable.h"
#include <QMessageBox>
#include <QtGlobal>

ProRataPeptideTable::ProRataPeptideTable( ProteomeInfo * ptmInfoInst, QWidget * qwParent ) 
	: ProRataTable( ptmInfoInst, qwParent )
{
	iRowCount = 0;
	iCurrentColumnSorted = 0;
	bIsCurrentSortAscending = true;
	connect( qtwTable, SIGNAL( clicked( const QModelIndex & ) ), 
			this, SLOT( rowClicked( const QModelIndex & ) ) );
	connect( qtwTable->horizontalHeader(), SIGNAL( sectionClicked( int ) ),
		this, SLOT( colClickedForSorting( int ) ) );

	qtwTable->horizontalHeader()->setSortIndicatorShown( true );
	qtwTable->horizontalHeader()->setSortIndicator( iCurrentColumnSorted, Qt::AscendingOrder );
	qtwTable->horizontalHeader()->resizeSection( 0, 180 );
	//qtwTable->horizontalHeader()->setSortIndicator( 0, Qt::AscendingOrder );
	
	ppepRatio = NULL;
}

ProRataPeptideTable::~ProRataPeptideTable()
{

}

void ProRataPeptideTable::setXmlProcessor( ProRataXmlProcessors * prxpProc )
{


	prxpProcessor = prxpProc;
	prxpProcessor->resetPointers();

	

	//this->populateTable();



}

void ProRataPeptideTable::rowClicked( const QModelIndex & qmiClickedItem )
{
	int iRow = qmiClickedItem.row();

	//QMessageBox::information( this, "selected row by mouse", QString::number( iRow ) );

	if ( iRow < 0 )
	{
		return;
	}

	qtwTable->selectRow( iRow ) ;

	delete ppepRatio;
	ppepRatio = new PeptideRatio;
	vpepInfo.at(iRow)->setPeptideRatio( ppepRatio );

	emit peptideClicked( ppepRatio );
	emit peptideClicked( vpepInfo.at(iRow) );
}
/*
void ProRataPeptideTable::rectMLESelected( const QwtDoubleRect rect )
{
	qDebug( " inside MLE " );
	int iRow = 0;
	bool bPeptideSelected = false;
	for( unsigned int i = 0; i < vpepInfo.size(); i++ )
	{
		if( !vpepInfo.at(i)->getValidity() )
			continue;
		if( rect.contains( vpepInfo.at(i)->getPCALog2Ratio(), vpepInfo.at(i)->getPCALog2SNR() ) )
		{
			bPeptideSelected = true;
			iRow = i;
			break;
		}
	}

	if( !bPeptideSelected )
		return;
	
//	QMessageBox::information( this, "Row to be selected", QString::number( iRow ) );
	qtwTable->selectRow( iRow ) ;
	//qtwTable->showRow( iRow );

	delete ppepRatio;
	ppepRatio = new PeptideRatio;
	vpepInfo.at(iRow)->setPeptideRatio( ppepRatio );

	emit peptideClicked( ppepRatio );
	emit peptideClicked( vpepInfo.at(iRow) );
}
*/

void ProRataPeptideTable::pointMLESelected( const QwtDoublePoint point )
{
	int iRow = 0;
	unsigned int iMinIndex = 0;
	double dMinDistance = 1000;
	double dCurrentDistance = 0;
	map< unsigned int, double > mIndexDistance;
	for( unsigned int i = 0; i < vpepInfo.size(); i++ )
	{
		if( !vpepInfo.at(i)->getValidity() )
			continue;
		// the x and y axes are not in the same range
		// take a sqrt to compress the difference
		dCurrentDistance = sqrt( fabs( point.x() - vpepInfo.at(i)->getPCALog2Ratio() ) )
			+ sqrt( fabs( point.y() - vpepInfo.at(i)->getPCALog2SNR() ) );
		if( dCurrentDistance < dMinDistance )
		{
			dMinDistance = dCurrentDistance;
			iMinIndex = i;
		}
	}

	if( dMinDistance > 1000 )
		return;

	//QMessageBox::information( this, "Row to be selected, before intising", QString::number( iMinIndex ) );
	iRow = (int)iMinIndex;
	
	//QMessageBox::information( this, "Row to be selected", QString::number( iRow ) );
	//qtwTable->selectRow( iRow ) ;
	//qtwTable->showRow( iRow );

    QAbstractItemView::SelectionMode sm = qtwTable->selectionMode();
    qtwTable->setSelectionMode(QAbstractItemView::ExtendedSelection);
    qtwTable->selectRow( iRow );
    qtwTable->setSelectionMode(sm);

	delete ppepRatio;
	ppepRatio = new PeptideRatio;
	vpepInfo.at(iRow)->setPeptideRatio( ppepRatio );

	emit peptideClicked( ppepRatio );
	emit peptideClicked( vpepInfo.at(iRow) );

	//qtwTable->selectRow( 0 ) ;
}

void ProRataPeptideTable::newProteinClicked( ProteinInfo * pproInfo )
{
	if (!(pproInfo))
	{
		vpepInfo.clear();
		populateTable();
		return;
	}
	vpepInfo.clear();
	vpepInfo = pproInfo->getPeptideInfo();	
	//populateTable();
//	qtwTable->horizontalHeader()->setSortIndicator( iCurrentColumnSorted, Qt::AscendingOrder );
	//qtwTable->horizontalHeader()->setSortIndicator( 0, Qt::AscendingOrder );
	qtwTable->horizontalHeader()->setStretchLastSection( true );
//	colClickedForSorting( iCurrentColumnSorted );
	sortColumn( iCurrentColumnSorted, bIsCurrentSortAscending );


}

void ProRataPeptideTable::populateTable()
{

	qaiModel->removeRows( 0, qaiModel->rowCount( QModelIndex() ), QModelIndex() );
	iRowCount = 0;

	int iTupleSize = vpepInfo.size();

//	qaiModel->beginInsertRows();
	qaiModel->insertRows( iRowCount, iTupleSize, QModelIndex() );
//	qaiModel->endInsertRows();

	
	for( int i = 0; i < iTupleSize; i++ )
	{
		// Add Sequence
		qaiModel->setData( qaiModel->index( iRowCount, 0, QModelIndex() ),
				vpepInfo.at(i)->getSequence().c_str() );

		// Add Log2Ratio
		qaiModel->setData( qaiModel->index( iRowCount, 1, QModelIndex() ),
				QString::number( vpepInfo.at(i)->getPCALog2Ratio()) );

		// Add Eigenvalue Ratio
		qaiModel->setData( qaiModel->index( iRowCount, 2, QModelIndex() ),
				QString::number(  vpepInfo.at(i)->getPCALog2SNR() ) );

		// Add Validity
		qaiModel->setData( qaiModel->index( iRowCount, 3, QModelIndex() ),
				vpepInfo.at(i)->getValidity() );

		qtwTable->resizeRowToContents( iRowCount );
		iRowCount++;
	}
	//qtwTable->sortItems( 0 );
	//qaiModel->sort( 0 );
}



void ProRataPeptideTable::colClickedForSorting( int iColumn )
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

void ProRataPeptideTable::sortColumn( int iColumn, bool bIsAscending )
{
	if( bIsAscending )
	{
		iCurrentColumnSorted = iColumn;
		bIsCurrentSortAscending = true;
		qtwTable->horizontalHeader()->setSortIndicator( iColumn, Qt::AscendingOrder );
		switch ( iColumn )
		{
			case (0):
				mainProteomeInfoInstance->sortPeptideInfo( vpepInfo, "sequence" );
				break;
			case (1):
				mainProteomeInfoInstance->sortPeptideInfo( vpepInfo, "log2Ratio" );
				break;
			case (2):
				mainProteomeInfoInstance->sortPeptideInfo( vpepInfo, "log2SN" );
				break;
			case (3):
				mainProteomeInfoInstance->sortPeptideInfo( vpepInfo, "validity" );
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
			mainProteomeInfoInstance->sortPeptideInfoDescending( vpepInfo, "sequence" );
			break;
		case (1):
			mainProteomeInfoInstance->sortPeptideInfoDescending( vpepInfo, "log2Ratio" );
			break;
		case (2):
			mainProteomeInfoInstance->sortPeptideInfoDescending( vpepInfo, "log2SN" );
			break;
		case (3):
			mainProteomeInfoInstance->sortPeptideInfoDescending( vpepInfo, "validity" );
			break;
		default:
			break;
		}

	}
	emit flushGraph();
	populateTable();

}

void ProRataPeptideTable::cleanUp()
{
	vpepInfo.clear();
	populateTable();
	emit flushGraph();
	return;

}

void ProRataPeptideTable::keyPressEvent ( QKeyEvent * event )
{
		//QMessageBox::information( this, "selected row by keyboard", QString::number( qtwTable->currentIndex().row() ) );
}
