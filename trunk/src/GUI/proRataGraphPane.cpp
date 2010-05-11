
#include "proRataGraphPane.h"

#include <QPainter>
#include <QPrintDialog>
#include <QPrinter>
#include <QPixmap>
#include <QMessageBox>
#include <QFileDialog>
#include "helpIcons.h"

ProRataGraphPane::ProRataGraphPane( QWidget * qwParent )
		: QWidget( qwParent )
{

	/*
	pripProteinPane = new ProRataImagePane;

	priLikehoodImage = new ProRataImage( "../images/Radis1.jpeg" );
	priLikehoodImage->setTitle( "Likehood Curve" );

	priCoverageImage = new ProRataImage( "../images/Radis2.jpeg" );
	priCoverageImage->setTitle( "Sequence Coverage" );

	pripProteinPane->addImage( priLikehoodImage );
	pripProteinPane->addImage( priCoverageImage );

	pripPeptidePane = new ProRataImagePane;

	priChromatogram = new ProRataImage( "../images/Radis3.jpeg" );
	priChromatogram->setTitle( "Chromatogram" );

	priPCA = new ProRataImage( "../images/Radis4.jpeg" );
	priPCA->setTitle( "Prinicipal Component Analysis" );

	pripPeptidePane->addImage( priChromatogram );
	pripPeptidePane->addImage( priPCA );

	qsPaneDivider = new QSplitter;
	qsPaneDivider->setOrientation( Qt::Vertical );
	
	qsPaneDivider->addWidget( pripProteinPane );
	qsPaneDivider->addWidget( pripPeptidePane );

	qvbMainLayout = new QVBoxLayout;
	qvbMainLayout->addWidget( qsPaneDivider );
	*/

	buildMenu();

	qvbMainLayout = new QVBoxLayout;

	qtwGraphTab = new QTabWidget;
//	connect( qtwGraphTab, SIGNAL( currentChanged( int ) ),
//		this, SLOT( currentSelected(int) ) );
	qvbMainLayout->addWidget( qtwGraphTab );

	qszDialogSize.setWidth( 600 );
	qszDialogSize.setHeight( 600 );

	qpDialogPosition.setX( 0 );
	qpDialogPosition.setY( 0 );

//	qiSeen = QIcon( QPixmap( nonUpdated_xpm ) );
//	qiUpdated = QIcon( QPixmap( updated_xpm ) );

	/*
	qszDialogSize = new QSize( 600, 600 );
	qpDialogPosition = new QPoint( 0, 0 );
	*/

	setLayout( qvbMainLayout );
}

ProRataGraphPane::~ProRataGraphPane()
{

}

void ProRataGraphPane::addGraph( const QString qsName, ProRataGraph * prgGraph )
{
	/*
	if ( qsName == "Chromatogram Graph" )
	{
		cout << "Yess Chro got" << endl;
		connect( this, SIGNAL( selected(const QwtDoublePoint &pos) ), 
				prgGraph, SLOT( leftPressed(const QwtDoublePoint &pos) ) );
	}
	*/
	prgGraph->setPosition( qtwGraphTab->count() );
	vprgGraphs.push_back( prgGraph );
	//qtwGraphTab->addTab( prgGraph, qiSeen, qsName );
//	qtwGraphTab->insertTab( qtwGraphTab->count(), prgGraph, qiSeen, qsName );
	qtwGraphTab->insertTab( qtwGraphTab->count(), prgGraph, qsName );
//	connect( prgGraph, SIGNAL( updated(QString,int) ), this, SLOT( graphUpdated(QString,int) ) );

}

void ProRataGraphPane::contextMenuEvent( QContextMenuEvent * event )
{
	//cout << "Inside context menu" << endl;
	event->accept();
	qmContextMenu->exec( event->globalPos() );
}

void ProRataGraphPane::buildMenu()
{

	qmContextMenu = new QMenu( this );

	qaDetach = new QAction( "Detach", this );
	connect( qaDetach, SIGNAL( triggered() ), this, SLOT( detach() ) );

	// the print function is commented out.
	// Removing the comment signs will add the print function back in.
	
	
	qaExport = new QAction( "Export...", this );
	connect( qaExport, SIGNAL( triggered() ), this,
		SLOT( exportGraph() ) );
	

	qmContextMenu->addAction( qaDetach );
//	qmContextMenu->addSeparator();
	qmContextMenu->addAction( qaExport );
	
}

void ProRataGraphPane::detach()
{



	if ( !(qtwGraphTab->currentWidget()) )
	{
		return ;
	}

	//QString qsTitle( "ProRata: " );
	QString qsTitle("");
	int iCurrentIndex = qtwGraphTab->currentIndex();
	QString qsName;

	if( iCurrentIndex >= 0 )
	{
		qsName = qtwGraphTab->tabText( iCurrentIndex );
		qsTitle = qsTitle + qsName;
	}
	//cout << "Name = " << qsTitle.toAscii().data() << endl;

	ProRataGraph *clickedGraph = (ProRataGraph *)qtwGraphTab->currentWidget();
	qtwGraphTab->removeTab( iCurrentIndex );

//	ProRataGraphDialog *prgdChildDialog = new ProRataGraphDialog( this, Qt::Dialog | Qt::WindowStaysOnTopHint );
	ProRataGraphDialog *prgdChildDialog = new ProRataGraphDialog( this, Qt::Dialog );

	connect( prgdChildDialog, SIGNAL( closeAction( ProRataGraphDialog * ) ), 
			this, SLOT( closeDialog( ProRataGraphDialog * ) ) );

	connect( prgdChildDialog, SIGNAL( resizeSignal( const QSize & ) ), 
			this, SLOT( resizeSlot( const QSize & ) ) );

	connect( prgdChildDialog, SIGNAL( moveSignal( const QPoint & ) ), 
			this, SLOT( moveSlot( const QPoint & ) ) );

	//prgdChildDialog->setGraph( qsTitle, (ProRataGraph *)qtwGraphTab->currentWidget() );
	prgdChildDialog->setGraph( qsTitle, clickedGraph );
	prgdChildDialog->move( qpDialogPosition );
	prgdChildDialog->resize( qszDialogSize );

	/*
	Qt::WindowFlags flags = 0;
	flags = Qt::Dialog;
	flags |= Qt::WindowStaysOnTopHint;
	prgdChildDialog->setWindowFlags(flags);


	QPoint pos = prgdChildDialog->pos();
	if (pos.x() < 0)
		pos.setX(0);
	if (pos.y() < 0)
		pos.setY(0);
	prgdChildDialog->move(pos);
	*/
	prgdChildDialog->show();



	/*
	   QString qsTitle( "ProRata: " );
	   int iCurrentIndex = qtwGraphTab->currentIndex();
	QString qsName;

	if( iCurrentIndex >= 0 )
	{
		qsName = qtwGraphTab->tabText( iCurrentIndex );
		qsTitle = qsTitle + qsName;
	}
	cout << "Name = " << qsTitle.toAscii().data() << endl;

	if ( qtwGraphTab->currentWidget() )
	{
		QWidget *qwNewChildWindow = new QWidget( NULL, Qt::Dialog );
		qwNewChildWindow->setWindowTitle( qsTitle );
		QHBoxLayout *qhblChildMainLayout = new QHBoxLayout;
		//ProRataGraph *theGraph = getGraph( qsName );
		ProRataGraph *theGraph = new ProRataGraph;
		theGraph = getGraph( qsName );
		if ( theGraph )
		{
			//qhblChildMainLayout->addWidget( qtwGraphTab->currentWidget() );
			qhblChildMainLayout->addWidget( theGraph );
		}
		qwNewChildWindow->setLayout( qhblChildMainLayout );
		qwNewChildWindow->show();
	}
	*/

	/*
	QDialog * qdGraphHolder = new QDialog( this );
	qdGraphHolder->setModal( false );
	qdGraphHolder->setExtension( qtwGraphTab->currentWidget() );
	qdGraphHolder->show();
	*/
}

/*

ProRataGraph * ProRataGraphPane::getGraph( const QString &qsName )
{
	vector< ProRataGraph *>::iterator vprgpIter;

	for ( vprgpIter = vprgGraphs.begin();
			vprgpIter != vprgGraphs.end();
			vprgpIter++ )
	{
		if ( qsName == (*vprgpIter)->accessibleName() )
			cout << "Got it" << endl;
		return (*vprgpIter);
	}
	return NULL;
}
*/

void ProRataGraphPane::closeDialog( ProRataGraphDialog * graph )
{
	//cout << "Inside Dialog slot" << endl;
	int iPos = graph->getGraph()->getPosition();
	for (int i = 0; i < qtwGraphTab->count(); i++)
	{
		if ( iPos < ((ProRataGraph *)(qtwGraphTab->widget(i)))->getPosition() )
		{
		//	qtwGraphTab->insertTab( i, graph->getGraph(), qiSeen, graph->getGraphName() );	
			qtwGraphTab->insertTab( i, graph->getGraph(), graph->getGraphName() );	
			return;
		}
	}
//	qtwGraphTab->insertTab( iPos, graph->getGraph(), qiSeen, graph->getGraphName() );	
	qtwGraphTab->insertTab( iPos, graph->getGraph(), graph->getGraphName() );	
}

void ProRataGraphPane::resizeSlot( const QSize &qszSz )
{
	//cout << "Setting new size" << endl;
	qszDialogSize = qszSz;	
}

void ProRataGraphPane::moveSlot( const QPoint &qptPt )
{
	//cout << "Setting new position" << endl;
	qpDialogPosition = qptPt;
}

/*
void ProRataGraphPane::newDataAdded()
{
}
*/

void ProRataGraphPane::print()
{
	if ( !(qtwGraphTab->currentWidget()) )
	{
		return ;
	}

/*
	QString qsSaveFileName = QFileDialog::getSaveFileName(
                    this,
                    "ProRata - Choose a file name to save",
                    ".",
                    "Images (*.png *.jpg)");

	if (qsSaveFileName.isEmpty())
	{
		return;
	}
	*/

	ProRataGraph *clickedGraph = (ProRataGraph *)qtwGraphTab->currentWidget();

	QPrinter printer;
    //printer.setOutputToFile(true);
    //printer.setOutputFileName(qsSaveFileName);
    printer.setDocName(QString("ProRata_Graph"));

    printer.setCreator("ProRata Graph");
    printer.setOrientation(QPrinter::Landscape);

	QPrintDialog dialog(&printer);
    if ( dialog.exec() )
    {
		QwtPlotPrintFilter filter;
		if ( printer.colorMode() == QPrinter::GrayScale )
        {
            filter.setOptions(QwtPlotPrintFilter::PrintAll 
                & ~QwtPlotPrintFilter::PrintCanvasBackground);
        }
		clickedGraph->print( printer, filter );
	}
	
	

/*
	//QPixmap *qpPic = new QPixmap( qsSaveFileName );
	//QPainter * qpntrPainter = new QPainter( qpPic );
	QPixmap qpPic( qsSaveFileName );
	qpPic.save( qsSaveFileName );
	//QPainter * qpntrPainter = new QPainter( qpPic );

	//clickedGraph->print( qpntrPainter, clickedGraph->frameGeometry() );
	clickedGraph->print( qpPic );
*/
	//QMessageBox::information( this, "Printing done.", QString( "Saved \"") + qsSaveFileName + QString( "\"" ) );

}

void ProRataGraphPane::exportGraph()
{
	ProRataGraph *clickedGraph = (ProRataGraph *)qtwGraphTab->currentWidget();

	if ( !( clickedGraph ) )
	{
		return ;
	}

	
	QString qsFName = QFileDialog::getSaveFileName( this,
			tr( "Choose a file to export the graph."), 
			QDir::homePath(), "PNG Files (*.png)" );
	if ( qsFName.isEmpty() )
	{       
		return;
	}

	if (!(qsFName.endsWith( QString( ".png" ), Qt::CaseInsensitive ) ) )
	{
		qsFName = qsFName + QString( ".png" );
	}

	clickedGraph->saveToFile( qsFName, "PNG" );

}

/*
void ProRataGraphPane::graphUpdated(QString qsName, int iCleanEvent)
{
	for (int i = 0; i < qtwGraphTab->count() ; i++)
	{
		if ( qtwGraphTab->tabText(i) == qsName)
		{
			//qtwGraphTab->setTabText( i, QString("<B>") + qsName + QString("</B>") );
			if (i == qtwGraphTab->currentIndex())
			{
				qtwGraphTab->setTabIcon( i, qiSeen );
				return;
			}
			qtwGraphTab->setTabIcon( i, 
				(iCleanEvent ? qiUpdated : qiSeen ) );

		}
	}
	
}

void ProRataGraphPane::currentSelected(int iIndex)
{
	qtwGraphTab->setTabIcon( iIndex, qiSeen );
}
*/
