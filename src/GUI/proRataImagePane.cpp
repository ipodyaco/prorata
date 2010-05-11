
#include "proRataImagePane.h"

ProRataImagePane::ProRataImagePane( QWidget * qwPane ) : QWidget( qwPane )
{
	qvbMainLayout = new QHBoxLayout;
	qvbMainLayout->setMargin( 0 );
	qvbMainLayout->setSpacing( 0 );

	qsImageDivider = new QSplitter;
	qsImageDivider->setOrientation( Qt::Horizontal );

	qvbMainLayout->addWidget( qsImageDivider );

	setLayout( qvbMainLayout );
}

ProRataImagePane::~ProRataImagePane()
{

}

void ProRataImagePane::addImage( ProRataImage * priNewGraph )
{
	vpriGraphs.push_back( priNewGraph );
	qsImageDivider->addWidget( priNewGraph );
	qsImageDivider->setStretchFactor( ( vpriGraphs.size() - 1 ), 1 );
}
