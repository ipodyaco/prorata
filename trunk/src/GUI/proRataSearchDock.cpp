
#include "proRataSearchDock.h"


ProRataSearchDock::ProRataSearchDock( QWidget* parent, Qt::WFlags fl )
    : QDockWidget( parent, fl )
{
	buildUI();
}

ProRataSearchDock::~ProRataSearchDock()
{
}

void ProRataSearchDock::buildUI()
{

	QWidget *qwTopLevel = new QWidget(this);

	findLayout = new QHBoxLayout;
	findLayout->setMargin( 0 );
	findLayout->setSpacing( 0 );

	qlFindLabel = new QLabel;
	qlFindLabel->setText( tr("Find:" ) );
	findLayout->addWidget( qlFindLabel );

	qleSearchString = new QLineEdit;
	qleSearchString->setMinimumWidth( 40 );
	findLayout->addWidget( qleSearchString );

	qpbSearchButton = new QPushButton;
	qpbSearchButton->setText( tr( "Search" ) );
	findLayout->addWidget( qpbSearchButton );

	qpbShowAll = new QPushButton;
	qpbShowAll->setText( tr( "Show All" ) );
	findLayout->addWidget( qpbShowAll );
	findLayout->addStretch();

	qwTopLevel->setLayout( findLayout );
	setWidget( qwTopLevel );
}
