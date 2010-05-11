
#include "proRataTablePane.h"
#include <QList>

ProRataTablePane::ProRataTablePane( ProRataXmlProcessors * prxpProcessor, 
		QWidget *qwPane ) : QWidget( qwPane )
{
	buildUI();
}

ProRataTablePane::ProRataTablePane( QWidget *qwPane ) : QWidget( qwPane )
{
	/*
	vector< ProteinInfo * > vpProteinInfo = mainProteomeInfo->getProteinInfo( "" );

	//Create a table instance.
	prtProtein = new ProRataProteinTable( mainProteomeInfo );

	// Create protein table title.
	prtProtein->setTableTitle( QString( "Protein Quantification Results" ) );

	// Create header info for the protein table.
	QStringList qslProteinTableHeader;
	qslProteinTableHeader << "Locus" <<
		"Log2 Ratio" << "Confidence Interval Width" <<
		"Description";

	// Setup the model by passing the table header information.
	prtProtein->setupModel( qslProteinTableHeader );

	//prtProtein->setXmlProcessor( prxpProcessor );
	prtProtein->populateTable( vpProteinInfo );
	*/

	/*
	prtPeptide = new ProRataPeptideTable( mainProteomeInfo );
	prtPeptide->setTableTitle( QString( "Peptide Quantification Results" ) );
	QStringList qslPeptideTableHeader;
	qslPeptideTableHeader << "Sequence" <<
		"Log2 Ratio" << "Eigenvalue Ratio" <<
		"Validity";
	prtPeptide->setupModel( qslPeptideTableHeader );
	*/
	//prtPeptide->setXmlProcessor( prxpProcessor );

	/*
	connect( prtProtein, SIGNAL( proteinClicked( const QString & ) ),
			prtPeptide, SLOT( newProteinClicked( const QString & ) ) );
			*/

	//prspFind = new ProRataSearchPane;
	//connect( prspFind, SIGNAL( searchStringSignal( QString ) ),
	//	this, SLOT( reEmitSlot( QString) ) );

	qvbLayout = new QVBoxLayout;
	qvbLayout->setSpacing( 0 );
	qvbLayout->setMargin( 0 );
	
//	QSplitter *qsTableDivider;
	qsTableDivider = new QSplitter( this );
	qsTableDivider->setOrientation( Qt::Vertical );
	//qvbLayout->addWidget( prspFind );
	//qsTableDivider->addWidget( prspFind );

	qvbLayout->addWidget( qsTableDivider );
	
	//qvbLayout->addWidget( prtProtein );
	//qvbLayout->addWidget( prtPeptide );

	setLayout( qvbLayout );
}

ProRataTablePane::~ProRataTablePane()
{

}

void ProRataTablePane::buildUI()
{
#if 0

	//Stuff for Protein Table.
	//Create a table instance.
	prtProtein = new ProRataProteinTable();

	// Create protein table title.
	prtProtein->setTableTitle( QString( "Protein Quantification Results" ) );

	// Create header info for the protein table.
	QStringList qslProteinTableHeader;
	qslProteinTableHeader << "Locus" <<
		"Log2 Ratio" << "Confidence Interval Width" <<
		"Description";

	// Setup the model by passing the table header information.
	prtProtein->setupModel( qslProteinTableHeader );

	//prtProtein->setXmlProcessor( prxpProcessor );

	prtPeptide = new ProRataPeptideTable();
	prtPeptide->setTableTitle( QString( "Peptide Quantification Results" ) );
	QStringList qslPeptideTableHeader;
	qslPeptideTableHeader << "Sequence" <<
		"Log2 Ratio" << "Eigenvalue Ratio" <<
		"Validity";
	prtPeptide->setupModel( qslPeptideTableHeader );
	//prtPeptide->setXmlProcessor( prxpProcessor );

	connect( prtProtein, SIGNAL( proteinClicked( const QString & ) ),
			prtPeptide, SLOT( newProteinClicked( const QString & ) ) );

	//prspFind = new ProRataSearchPane;

	qvbLayout = new QVBoxLayout;
	qvbLayout->setSpacing( 0 );
	qvbLayout->setMargin( 0 );
	
	/*
	qsTableDivider = new QSplitter( this );
	qsTableDivider->setOrientation( Qt::Vertical );
	qsTableDivider->addWidget( prtProtein );
	qsTableDivider->addWidget( prtPeptide );
	qsTableDivider->setStretchFactor( 0, 1 );
	qsTableDivider->setStretchFactor( 1, 1 );


	qvbLayout->addWidget( qsTableDivider );
	*/

	qvbLayout->addWidget( prspFind );
	qvbLayout->addWidget( prtProtein );
	qvbLayout->addWidget( prtPeptide );


	setLayout( qvbLayout );
#endif
}

void ProRataTablePane::setProteinTable( ProRataTable * prtTable )
{	
	prtProtein = (ProRataProteinTable*) prtTable;	
	qsTableDivider->addWidget( prtProtein );
	
	QList<int> qlSizes;
	//qlSizes.push_back( 15 );
	qlSizes.push_back( 475 );
	qlSizes.push_back( 300 );

	qsTableDivider->setSizes( qlSizes );
}


void ProRataTablePane::setPeptideTable( ProRataTable * prtTable )
{	
	prtPeptide = (ProRataPeptideTable*) prtTable;
	qsTableDivider->addWidget( prtPeptide );
}

/*
void ProRataTablePane::hideSearchPane( bool bStatus )
{
	if ( bStatus )
	{
		prspFind->hide();
	}
	else
	{
		prspFind->show();
	}
}

void ProRataTablePane::reEmitSlot( QString qsStr )
{
	qDebug( "Inside ProRataTablePane reEmit" );
	emit( qsStr );
}
*/
