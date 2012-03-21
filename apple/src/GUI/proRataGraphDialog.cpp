
#include "proRataGraphDialog.h"

ProRataGraphDialog::ProRataGraphDialog( QWidget * qwParent, Qt::WFlags flag )
		: QWidget( qwParent, flag )
{

	prgGraph = NULL;
	qsGraphName = "";
	qhblChildMainLayout = new QHBoxLayout;
	setLayout( qhblChildMainLayout );
}

ProRataGraphDialog::~ProRataGraphDialog()
{
}

void ProRataGraphDialog::setGraph( const QString qsName, ProRataGraph *graph )
{
	prgGraph = graph;
	qsGraphName = qsName;
	setWindowTitle( QString( "ProRata: " ) + qsName );
	qhblChildMainLayout->addWidget( graph );
}

void ProRataGraphDialog::closeEvent( QCloseEvent * /*event*/ )
{
	//cout << "Inside dialog's close event" << endl;
	emit closeAction( this );
}

void ProRataGraphDialog::resizeEvent( QResizeEvent * event )
{
	emit resizeSignal( event->size() );
}

void ProRataGraphDialog::moveEvent ( QMoveEvent * event )
{
	emit moveSignal( event->pos() );
}
