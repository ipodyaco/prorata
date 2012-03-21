
#include "proRataMainWindow.h"
#include <QMessageBox>
#include <QSplashScreen>
#include <stdio.h>
#include "proRataSplash.h"


ProRataMainWindow::ProRataMainWindow(QWidget *parent, Qt::WFlags flags)
        : QMainWindow(parent, flags)

{
	setWindowIcon( QIcon( QPixmap( proRataMainIcon ) ) );

	prMainWidget = NULL;

	qwCentralWidget = new QWidget;
	qvblCentralLayout = new QVBoxLayout;
	qvblCentralLayout->setSpacing( 0 );
	qvblCentralLayout->setMargin( 2 );

	createActions();
	createMenus();
	createToolBars();
	createStatusBar();

	readSettings();
	QDir::setCurrent( QDir::rootPath() );

	setWindowTitle( tr( "ProRata 1.1" ) );
	//paintEvent();

	showMaximized();
	//paintEvent();

	qwCentralWidget->setLayout( qvblCentralLayout );
	setCentralWidget( qwCentralWidget );
	//qwCentralWidget->hide();
	//QPixmap *qpImage = new QPixmap( proRataSplashScreen_xpm );

	/*
	QWidget *qwCredits = new QWidget( this, Qt::SplashScreen );
	QVBoxLayout *ly = new QVBoxLayout;
	qwCredits->setLayout( ly );
	ly->addWidget( new QIcon(*qpImage) );
	qwCredits->show();
	*/	

	//qssSplash =  new QSplashScreen( *qpImage, Qt::SplashScreen );
	//qssSplash->show();

}

void ProRataMainWindow::open()
{
	//delete qssSplash;
	int iAns = 2;
	if ( prMainWidget != NULL )
	{
		iAns = QMessageBox::question(
                this,
                tr("ProRata -- Closing current session."),
                tr("Do you want to close the current session?"),
                tr("&Yes"), tr("&No"),
                QString(), 0, 1);

		if ( iAns == 0)
		{
			mainClose();
		}
		else
		{
			return;
		}
	}

	QFileDialog *qfdResultsFile = new QFileDialog;

	QString qsFName = qfdResultsFile->getOpenFileName( this,
			"Choose a ProRata results file to open", NULL,
			"Results (*.qpr.xml)" );
	
	if ( qsFName.isEmpty() )
		return;
	
	statusBar()->showMessage( tr( "Loading..." ) );
	qpbProgress->setValue( 0 );

	qsFilename = qsFName;
	qfdResultsFile->selectFile( qsFName );

	QString qsDirectory = qsFName;
	int iIndex = qsDirectory.lastIndexOf( "/" );
	iIndex++;
	qsDirectory = qsDirectory.remove( iIndex, qsDirectory.length() - iIndex );

#if _WIN32
	qsDirectory.replace( "/", "\\" );
#endif

	QString qsConfigFileName = qsDirectory + "ProRataConfig.xml";
	if( !ProRataConfig::setFilename( qsConfigFileName.toAscii().data() ) )
	{
		QMessageBox::information(this, "Error", "Cannot load ProRataConfig.xml! Make sure ProRataConfig.xml is in the same directory as the qpr file" );
		return;
	}
	ProRataConfig::setWorkingDirectory( qsDirectory.toAscii().data() );

	prMainWidget = new ProRataMainWidget;
	connect( prMainWidget, SIGNAL( updateStatus( int ) ),
		qpbProgress, SLOT( setValue( int ) ) );
	connect( prMainWidget, SIGNAL( updateStatusMessage( QString ) ),
		statusBar(), SLOT( showMessage( QString ) ) );

	prMainWidget->setFilename( qsFilename );
	prMainWidget->redirectSearchSlot( qaSearch->isChecked() );
	qpbProgress->setValue( 10 );
	statusBar()->showMessage( tr( "Loading.." ) );
	qpbProgress->reset();
	statusBar()->showMessage( tr( "Done" ) );
	readSettings();

	qvblCentralLayout->addWidget( prMainWidget );
	qwCentralWidget->show();
}

void ProRataMainWindow::open(QString qsFName)
{

	if ( prMainWidget != NULL )
	{
		mainClose();
	}


	if ( qsFName.isEmpty() )
		return;	

	qsFilename = qsFName;

	QString qsDirectory = qsFName;
	int iIndex = qsDirectory.lastIndexOf( "/" );
	iIndex++;
	qsDirectory = qsDirectory.remove( iIndex, qsDirectory.length() - iIndex );

#if _WIN32
	qsDirectory.replace( "/", "\\" );
#endif

	QString qsConfigFileName = qsDirectory + "ProRataConfig.xml";
	if( !ProRataConfig::setFilename( qsConfigFileName.toAscii().data() ) )
	{
		QMessageBox::information(this, "Error", "Cannot load ProRataConfig.xml! Make sure ProRataConfig.xml is in the same directory as the qpr file" );
		return;
	}
	ProRataConfig::setWorkingDirectory( qsDirectory.toAscii().data() );

	prMainWidget = new ProRataMainWidget(qwCentralWidget);
	prMainWidget->setFilename( qsFilename );
	qvblCentralLayout->addWidget( prMainWidget );
	prMainWidget->redirectSearchSlot( qaSearch->isChecked() );

}

void ProRataMainWindow::print()
{
}

void ProRataMainWindow::mainClose()
{
	delete prMainWidget;
	prMainWidget = NULL;
	readSettings();

}

void ProRataMainWindow::quit()
{
}


void ProRataMainWindow::defaultView()
{

	if ( prMainWidget )
	{
		prMainWidget->showTables( true );
		prMainWidget->showText( true );
		prMainWidget->showGraphs( true );
		qaShowTablesAct->setChecked( true );
		qaShowTextAct->setChecked( true );
		qaShowGraphsAct->setChecked( true );
		prMainWidget->updateUI();
	}
}

void ProRataMainWindow::showTables()
{
	if ( prMainWidget )
	{
		if ( qaShowTablesAct->isChecked() )
		{
			prMainWidget->showTables( true );
		}
		else
		{
			prMainWidget->showTables( false );
		}

	}
}

void ProRataMainWindow::showText()
{
	if ( prMainWidget )
	{
		if ( qaShowTextAct->isChecked() )
		{
			prMainWidget->showText( true );
		}
		else
		{
			prMainWidget->showText( false );
		}

	}
}

void ProRataMainWindow::showGraphs()
{
	if ( prMainWidget )
	{
		if ( qaShowGraphsAct->isChecked() )
		{
			prMainWidget->showGraphs( true );
		}
		else
		{
			prMainWidget->showGraphs( false );
		}

	}
}

void ProRataMainWindow::executeProRata()
{

	ProRataExecDialog * predStart = new ProRataExecDialog(this );
	connect( predStart, SIGNAL( qprFilename( QString ) ), this, SLOT( open( QString ) ) );
	//predStart->setAttribute( Qt::WA_ShowModal );
	predStart->setModal( true );
	predStart->exec();

}

void ProRataMainWindow::search()
{
	if ( prMainWidget != NULL )
	{
		prMainWidget->redirectSearchSlot( qaSearch->isChecked() );
	}
}

void ProRataMainWindow::configuration()
{
	ProRataConfigDialog * prcdConfig = new ProRataConfigDialog(this,Qt::Dialog);
	prcdConfig->setAttribute( Qt::WA_ShowModal );
	prcdConfig->show();
}

void ProRataMainWindow::mzXMLConversion()
{
	ProRataRaw2MzXMLBrowser * prr2mStart = new ProRataRaw2MzXMLBrowser(this);
	prr2mStart->setAttribute( Qt::WA_ShowModal );
	prr2mStart->show();

}

void ProRataMainWindow::mergeDirId()
{
	ProRataMerge * prmStart = new ProRataMerge(this );
	prmStart->setAttribute( Qt::WA_ShowModal );
	prmStart->exec();

}

void ProRataMainWindow::contents()
{
}

void ProRataMainWindow::generalHelp()
{
}

void ProRataMainWindow::introduction()
{
}

void ProRataMainWindow::about()
{
}

void ProRataMainWindow::createActions()
{
	qaOpenAct = new QAction( QIcon( QPixmap( open_xpm ) ), tr( "&Open..." ), this );
	qaOpenAct->setShortcut( tr( "Ctrl+O" ) );
	qaOpenAct->setStatusTip( tr ( "Open an existing ProRata Results file." ) );
	connect( qaOpenAct, SIGNAL( triggered() ), this, SLOT( open() ) );

	qaPrintAct = new QAction( QIcon( QPixmap( print_xpm ) ), tr( "&Print..." ), this );
	qaPrintAct->setStatusTip( tr( "Print the graphs."));
	connect(qaPrintAct, SIGNAL( triggered() ), this, SLOT( print() ) );
	qaCloseAct = new QAction( tr( "&Close" ), this );
	qaCloseAct->setStatusTip( tr( "Close the current Results file." ) );
	connect(qaCloseAct, SIGNAL( triggered() ), this, SLOT( mainClose() ) );

	qaQuitAct = new QAction( tr( "&Exit" ), this );
	qaQuitAct->setStatusTip( tr( "Quit ProRata" ) );
	connect(qaQuitAct, SIGNAL( triggered() ), this, SLOT( close() ) );

	qaDefaultViewAct = new QAction( tr( "&Default View" ), this );
	qaDefaultViewAct->setStatusTip( tr( "Revert to default view" ) );
	connect(qaDefaultViewAct, SIGNAL( triggered() ), this, SLOT( defaultView() ) );

	qaSearch = new QAction( QIcon( QPixmap( search_xpm ) ), tr( "&Search" ), this );
	qaSearch->setStatusTip( tr( "Search data" ) );
	qaSearch->setCheckable(true);
	qaSearch->setChecked(true);
	connect(qaSearch, SIGNAL( triggered() ), this, SLOT( search() ) );

	qaShowTablesAct = new QAction( tr( "Show &Tables" ), this );
	qaShowTablesAct->setStatusTip( tr( "Show/Hide the Protein and Peptide tables" ) );
	qaShowTablesAct->setCheckable( true );
	qaShowTablesAct->setChecked( true );
	connect(qaShowTablesAct, SIGNAL( triggered() ), this, SLOT( showTables() ) );

	qaShowTextAct = new QAction( tr( "Show Te&xt" ), this );
	qaShowTextAct->setStatusTip( tr( "Show/Hide textual information." ) );
	qaShowTextAct->setCheckable( true );
	qaShowTextAct->setChecked( true );
	connect(qaShowTextAct, SIGNAL( triggered() ), this, SLOT( showText() ) );

	qaShowGraphsAct = new QAction( tr( "Show &Graph" ), this );
	qaShowGraphsAct->setStatusTip( tr( "Show/Hide Graphs." ) );
	qaShowGraphsAct->setCheckable( true );
	qaShowGraphsAct->setChecked( true );
	connect(qaShowGraphsAct, SIGNAL( triggered() ), this, SLOT( showGraphs() ) );

	qaExecuteProRata = new QAction( tr( "&Execute ProRata..." ), this );
	qaExecuteProRata->setStatusTip( tr( "Start ProRata" ) );
	connect(qaExecuteProRata, SIGNAL( triggered() ), this, SLOT( executeProRata() ) );

	qaConfiguration = new QAction( tr( "&Create Config..." ), this );
	qaConfiguration->setStatusTip( tr( "Generate new ProRataConfig.xml." ) );
	connect(qaConfiguration, SIGNAL( triggered() ), this, SLOT( configuration() ) );

	qaMzXMLConvert = new QAction( tr( "&Raw 2 mzXML..." ), this );
	qaMzXMLConvert->setStatusTip( tr( "Convert RAW files to mzXML files." ) );
	connect(qaMzXMLConvert, SIGNAL( triggered() ), this, SLOT( mzXMLConversion() ) );

	qaMerge = new QAction( tr( "&Merge Directories..." ), this );
	qaMerge->setStatusTip( tr( "Merge Identification data directories." ) );
	connect(qaMerge, SIGNAL( triggered() ), this, SLOT( mergeDirId() ) );

	qaContentsAct = new QAction( tr( "&Contents" ), this );
	qaContentsAct->setStatusTip( tr( "Open main help pages." ) );
	connect(qaContentsAct, SIGNAL( triggered() ), this, SLOT( contents() ) );

	qaGeneralHelpAct = new QAction( tr( "&General Help" ), this );
	qaGeneralHelpAct->setStatusTip( tr( "Open help documents." ) );
	connect(qaGeneralHelpAct, SIGNAL( triggered() ), this, SLOT( generalHelp() ) );

	qaIntroductionAct = new QAction( tr( "&Introduction" ), this );
	qaIntroductionAct->setStatusTip( tr( "Start an introduction to ..." ) );
	connect(qaIntroductionAct, SIGNAL( triggered() ), this, SLOT( introduction() ) );

	qaAboutAct = new QAction( tr( "&About ProRata..." ), this );
	qaAboutAct->setStatusTip( tr( "About the software." ) );
	connect(qaAboutAct, SIGNAL( triggered() ), this, SLOT( about() ) );
}

void ProRataMainWindow::createMenus()
{
	qmFileMenu = menuBar()->addMenu( tr( "&File" ) );
	qmFileMenu->addAction( qaOpenAct );
	qmFileMenu->addSeparator();
	qmFileMenu->addAction( qaCloseAct );
	qmFileMenu->addAction( qaQuitAct );


	qmViewMenu = menuBar()->addMenu( tr( "&View" ) );
	qmViewMenu->addAction( qaDefaultViewAct );
	qmViewMenu->addSeparator();
	qmViewMenu->addAction( qaShowTablesAct );
	qmViewMenu->addAction( qaShowTextAct );
	qmViewMenu->addAction( qaShowGraphsAct );

	qmToolsMenu = menuBar()->addMenu( tr( "&Tools" ) );
	qmToolsMenu->addAction( qaExecuteProRata );
	qmToolsMenu->addSeparator();
	qmToolsMenu->addAction( qaConfiguration );
	qmToolsMenu->addAction( qaMzXMLConvert );
	qmToolsMenu->addAction( qaMerge );
	qmToolsMenu->addSeparator();
	qmToolsMenu->addAction( qaSearch );
	



	menuBar()->addSeparator();

	qmHelpMenu = menuBar()->addMenu( tr( "&Help" ) );
	qmHelpMenu->addAction( qaContentsAct );
	qmHelpMenu->addSeparator();
	qmHelpMenu->addAction( qaGeneralHelpAct );
	qmHelpMenu->addAction( qaIntroductionAct );
	qmHelpMenu->addSeparator();
	qmHelpMenu->addAction( qaAboutAct );
	qaContentsAct->setEnabled( false );
	qaGeneralHelpAct->setEnabled( false );
	qaIntroductionAct->setEnabled( false );
	qaAboutAct->setEnabled( false );

}

void ProRataMainWindow::createToolBars()
{
	qtbFrequentlyUsedActToolBar = addToolBar( tr( "FrequentlyUsed" ) );
	qtbFrequentlyUsedActToolBar->addAction( qaOpenAct );
	qtbFrequentlyUsedActToolBar->addSeparator();
	qtbFrequentlyUsedActToolBar->addAction( qaSearch );

}

void ProRataMainWindow::createStatusBar()
{
	qpbProgress = new QProgressBar;
	qpbProgress->setMaximumWidth( 250 );
	qpbProgress->setMaximumHeight( 15 );
	qpbProgress->setRange( 1, 10 );
	statusBar()->addPermanentWidget( qpbProgress );
	statusBar()->showMessage( tr( "Done" ) );
}

void ProRataMainWindow::readSettings()
{
	if (!prMainWidget)
	{
		qaCloseAct->setEnabled( false );
		qaShowTablesAct->setEnabled( false );
		qaShowTextAct->setEnabled( false );
		qaShowGraphsAct->setEnabled( false );
		qaSearch->setEnabled( false );
		qaDefaultViewAct->setEnabled( false );
	}
	else
	{

		qaCloseAct->setEnabled( true );
		qaShowTablesAct->setEnabled( true );
		qaShowTextAct->setEnabled( true );
		qaShowGraphsAct->setEnabled( true );
		qaSearch->setEnabled( true );
		qaDefaultViewAct->setEnabled( true );
	}
}

void ProRataMainWindow::loadFile( const QString & /*fileName */ )
{
}

void ProRataMainWindow::setCurrentFile( const QString & /*fileName*/ )
{
}

QString ProRataMainWindow::strippedName( const QString & /*fullFileName*/ )
{
	return QString( "" );
}

