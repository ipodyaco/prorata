
#include "proRataMainWidget.h"
#include <QList>
#include <QMessageBox>

ProRataMainWidget::ProRataMainWidget( const QString qsFname, QWidget * qwParent ) :
	QWidget( qwParent )
{
	emit updateStatus( 2 );
	emit updateStatusMessage( QString("Creating Tables") );

	prtpTables = NULL;
	prgpGraphs = NULL;
	prtpText = NULL;
	prspFind = NULL;
	mainProteomeInfo = NULL;
	mainMassSpecData = NULL;

	setFilename( qsFname );
}

ProRataMainWidget::ProRataMainWidget( QWidget * qwParent ) :
	QWidget( qwParent )
{
	prtpTables = NULL;
	prgpGraphs = NULL;
	prtpText = NULL;
	mainProteomeInfo = NULL;
	mainMassSpecData = NULL;
}

ProRataMainWidget::~ProRataMainWidget()
{
	delete mainProteomeInfo;
	delete mainMassSpecData;
}

void ProRataMainWidget::setFilename( const QString & qsFname )
{

	qsFilename = qsFname;

    mainProteomeInfo = new ProteomeInfo;
	mainProteomeInfo->readFileQPR( qsFname.toAscii().data() );

	mainMassSpecData = new ProRataMassSpectrumData;

	vector< ProteinInfo * > vpProteinInfo = mainProteomeInfo->getProteinInfo( "" );

	prtProtein = new ProRataProteinTable( mainProteomeInfo );
	prtProtein->setTableTitle( QString( "Protein Quantification Results" ) );
	QStringList qslProteinTableHeader;
	qslProteinTableHeader << "Locus"  << "Log2R" << "CIW" << "Description";
	prtProtein->setupModel( qslProteinTableHeader );
	prtProtein->populateTable( vpProteinInfo );

	prtPeptide = new ProRataPeptideTable( mainProteomeInfo );
	prtPeptide->setTableTitle( QString( "Peptide Quantification Results" ) );
	QStringList qslPeptideTableHeader;
	qslPeptideTableHeader << "Sequence" << "Log2R" << "Log2SNR" <<  "Validity";
	prtPeptide->setupModel( qslPeptideTableHeader );

	emit updateStatus( 2 );
	emit updateStatusMessage( QString("Populating Tables...") );

	connect( prtProtein, SIGNAL( flushGraph() ),
			prtPeptide, SLOT( cleanUp() ) );

	connect( prtProtein, SIGNAL( proteinClicked( ProteinInfo * ) ),
			prtPeptide, SLOT( newProteinClicked( ProteinInfo * ) ) );

	prtpTables = new ProRataTablePane;
	connect( prtpTables, SIGNAL( reEmitSearchString( QString ) ),
		this, SLOT( newSearchString( QString ) ) );

	prtpTables->setProteinTable( prtProtein );
	prtpTables->setPeptideTable( prtPeptide );

	emit updateStatus( 4 );


	// Creating Graphs.
	prgpGraphs = new ProRataGraphPane;

	ProRataLikelihoodGraph *likeliGraph = new ProRataLikelihoodGraph;
	connect( prtProtein, SIGNAL( flushGraph() ),
			likeliGraph, SLOT( cleanGraph() ) );
	connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
			likeliGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
	prgpGraphs->addGraph( QString( "Likelihood Plot"), likeliGraph );
	connect( likeliGraph, SIGNAL(  pointMLESelected( const QwtDoublePoint ) ),
			prtPeptide, SLOT( pointMLESelected( const QwtDoublePoint ) ) );
	connect( prtPeptide, SIGNAL( peptideClicked( PeptideInfo * ) ),
			likeliGraph, SLOT( peptideUpdated( PeptideInfo * ) ) );	

	if( QFile::exists(  QString::fromStdString( ProRataConfig::getFASTAfilename() ) ) )
	{
		// creat the sequence coverage plot
		ProRataSequenceCoverageGraph *coverageGraph = new ProRataSequenceCoverageGraph;
		connect( prtProtein, SIGNAL( proteinClicked( ProteinInfo * ) ),
				coverageGraph, SLOT( proteinUpdated( ProteinInfo * ) ) );
		connect( prtProtein, SIGNAL( flushGraph() ),
				coverageGraph, SLOT( cleanGraph() ) );
		connect( prtPeptide, SIGNAL( peptideClicked( PeptideInfo * ) ),
			coverageGraph, SLOT( peptideUpdated( PeptideInfo * ) ) );	
		prgpGraphs->addGraph( QString( "Seq Cov"), coverageGraph );
	}
	
	emit updateStatus( 5 );
	emit updateStatusMessage( QString("Building Graph Templates...") );


	QDir qDirMzXML( QString::fromStdString( ProRataConfig::getMZxmlDirectory() ) );
	QDir qDirXic( QString::fromStdString(ProRataConfig::getXICxmlDirectory() ) );
	
	if( qDirXic.exists() )
	{
		// creat ion chromatogram graph
		ProRataChromatogramGraph *chroGraph = new ProRataChromatogramGraph;
		connect( prtPeptide, SIGNAL( flushGraph() ),
				chroGraph, SLOT( cleanGraph() ) );
		connect( prtPeptide, SIGNAL( peptideClicked( PeptideRatio * ) ),
				chroGraph, SLOT( peptideUpdated( PeptideRatio * ) ) );
		connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
				chroGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
		prgpGraphs->addGraph( QString( "Ion Chromatogram"), chroGraph );

		// creat PPC chromatogram graph
		ProRataPPCGraph *ppcGraph = new ProRataPPCGraph;
		connect( prtPeptide, SIGNAL( flushGraph() ),
				ppcGraph, SLOT( cleanGraph() ) );
		connect( prtPeptide, SIGNAL( peptideClicked( PeptideRatio * ) ),
				ppcGraph, SLOT( peptideUpdated( PeptideRatio * ) ) );
		connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
				ppcGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
		prgpGraphs->addGraph( QString( "PPC Chromatogram"), ppcGraph );

		// creat PCA plot graph
		ProRataPCAGraph *pcaGraph = new ProRataPCAGraph;
		connect( prtPeptide, SIGNAL( flushGraph() ),
				pcaGraph, SLOT( cleanGraph() ) );
		connect( prtPeptide, SIGNAL( peptideClicked( PeptideRatio * ) ),
				pcaGraph, SLOT( peptideUpdated( PeptideRatio * ) ) );
		connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
				pcaGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
		connect( chroGraph, SIGNAL( MS1selected( long ) ),
				pcaGraph, SLOT( MS1selected( long ) ) );
		prgpGraphs->addGraph( QString( "PCA Plot"), pcaGraph );


		if( qDirMzXML.exists() )
		{
			// creat mass spectrum graph
			ProRataMassSpectrum *massSpecGraph = new ProRataMassSpectrum;
			massSpecGraph->setData( mainMassSpecData );
			connect( prtPeptide, SIGNAL( flushGraph() ),
					massSpecGraph, SLOT( cleanGraph() ) );
			connect( prtPeptide, SIGNAL( peptideClicked( PeptideRatio * ) ),
					massSpecGraph, SLOT( peptideUpdated( PeptideRatio * ) ) );
			connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
					massSpecGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
			connect( chroGraph, SIGNAL( MS1selected( long ) ),
					massSpecGraph, SLOT( updateMSGraph( long ) ) );
			prgpGraphs->addGraph( QString( "MS1 Scan"), massSpecGraph );

			// creat MS2 spectra graph
			ProRataTandemMassSpectrum *tandemMassSpecGraph = new ProRataTandemMassSpectrum;
			tandemMassSpecGraph->setData( mainMassSpecData );
			connect( prtPeptide, SIGNAL( flushGraph() ),
					tandemMassSpecGraph, SLOT( cleanGraph() ) );
			connect( prtPeptide, SIGNAL( peptideClicked( PeptideRatio * ) ),
					tandemMassSpecGraph, SLOT( peptideUpdated( PeptideRatio * ) ) );
			connect( prtProtein, SIGNAL( proteinClicked( ProteinRatio * ) ),
					tandemMassSpecGraph, SLOT( proteinUpdated( ProteinRatio * ) ) );
			connect( chroGraph, SIGNAL( MS2selected( long ) ),
					tandemMassSpecGraph, SLOT( updateMSGraph( long ) ) );
			prgpGraphs->addGraph( QString( "MS2 Scan"), tandemMassSpecGraph );
		}
	}

	emit updateStatus( 6 );

	// Creating Text pane.
	prtpText = new ProRataTextPane;
	prtpText->setProteinTable( prtProtein );
	prtpText->setPeptideTable( prtPeptide );
	
	QVBoxLayout *qvbLayoutMainLayout = new QVBoxLayout;

	QSplitter *qspGraphAndText = new QSplitter;
	qspGraphAndText->setOrientation( Qt::Vertical );
	qspGraphAndText->setOpaqueResize( false );

	QSplitter *qspMainSpliter = new QSplitter;
	qspMainSpliter->setOrientation( Qt::Horizontal );
	qspMainSpliter->setOpaqueResize( false );

	emit updateStatus( 7 );
	

	QList<int> qlSizes;
	qlSizes.push_back( 780 );
	qlSizes.push_back( 320 );
	
	qspGraphAndText->addWidget( prgpGraphs ); 
	qspGraphAndText->addWidget( prtpText );
	qspGraphAndText->setSizes( qlSizes );

	QList<int> qlSizes2;
	qlSizes2.push_back( 500 );
	qlSizes2.push_back( 800 );

	qwTableTop = new QWidget;
	QVBoxLayout * qvblTableTopLayout = new QVBoxLayout;

	qvblTableTopLayout->setMargin( 0 );
	qvblTableTopLayout->setSpacing( 0 );

	prspFind = new ProRataSearchPane;
	connect( prspFind, SIGNAL( searchStringSignal( QString ) ),
		this, SLOT( newSearchString( QString) ) );

	qvblTableTopLayout->addWidget( prspFind );
	qvblTableTopLayout->addWidget( prtpTables );

	qwTableTop->setLayout( qvblTableTopLayout );

	emit updateStatus( 8 );
	emit updateStatusMessage( QString("Rendering Workspace...") );


	qspMainSpliter->addWidget( qwTableTop );
	qspMainSpliter->addWidget( qspGraphAndText );
	qspMainSpliter->setSizes( qlSizes2 );

	qvbLayoutMainLayout->addWidget( qspMainSpliter );

	qvbLayoutMainLayout->setSpacing( 0 );
	qvbLayoutMainLayout->setMargin( 0 );

	setLayout( qvbLayoutMainLayout );
	emit updateStatus( 9 );

}

void ProRataMainWidget::showTables( bool bShow )
{
	if ( bShow )
	{
		//prtpTables->show();
		qwTableTop->show();
	}
	else
	{
		//prtpTables->hide();
		qwTableTop->hide();
	}
}

void ProRataMainWidget::showText( bool bShow )
{
	if ( bShow )
	{
		prtpText->show();
	}
	else
	{
		prtpText->hide();
	}
}

void ProRataMainWidget::showGraphs( bool bShow )
{
	if ( bShow )
	{
		prgpGraphs->show();
	}
	else
	{
		prgpGraphs->hide();
	}
}

void ProRataMainWidget::updateUI()
{
	prgpGraphs->update();

}

bool ProRataMainWidget::redirectSearchSlot( bool bStatus )
{
	if ( prspFind != NULL )
	{
		if ( bStatus )
		{
			prspFind->show();
		}
		else
		{
			prspFind->hide();
		}
	}

	return true;
}
		

void ProRataMainWidget::newSearchString( QString qsNewStr )
{
	//QMessageBox::information(this, "String string = ", qsNewStr.toAscii().data() );
	vector< ProteinInfo * > vpProteinInfo = mainProteomeInfo->getProteinInfo( qsNewStr.toAscii().data() );
	prtProtein->populateTable( vpProteinInfo );
}

