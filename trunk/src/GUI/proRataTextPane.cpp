
#include "proRataTextPane.h"
#include "proRataProteinTable.h"

ProRataTextPane::ProRataTextPane( QWidget * qwParent ) : 
	QWidget( qwParent )
{
	prtaText = new ProRataTextArea;
	qhblMainLayout = new QHBoxLayout;
	qhblMainLayout->addWidget( prtaText );
	setLayout( qhblMainLayout );
}

ProRataTextPane::~ProRataTextPane()
{

}
void ProRataTextPane::addTextArea( ProRataTextArea * prtaChild )
{
	
	qhblMainLayout->addWidget( prtaChild );
}

void ProRataTextPane::setProteinTable( ProRataProteinTable * prtTable )
{
	connect( prtTable, SIGNAL( proteinClicked( ProteinInfo * ) ),
		prtaText, SLOT( updateProteinDesc( ProteinInfo * ) ) );
	connect( prtTable, SIGNAL( proteinClicked( ProteinRatio * ) ),
		prtaText, SLOT( updateProteinDesc( ProteinRatio * ) ) );
	connect( prtTable, SIGNAL( flushGraph() ),
		prtaText, SLOT( cleanUp() ) );
}

void ProRataTextPane::setPeptideTable( ProRataPeptideTable * prtTable )
{
	connect( prtTable, SIGNAL( peptideClicked( PeptideInfo * ) ),
		prtaText, SLOT( updatePeptideDesc( PeptideInfo * ) ) );
	connect( prtTable, SIGNAL( peptideClicked( PeptideRatio * ) ),
		prtaText, SLOT( updatePeptideDesc( PeptideRatio * ) ) );
	connect( prtTable, SIGNAL( flushGraph() ),
		prtaText, SLOT( cleanUpPeptide() ) );
}

