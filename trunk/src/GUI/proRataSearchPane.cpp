
#include "proRataSearchPane.h"


ProRataSearchPane::ProRataSearchPane( QWidget* parent, Qt::WFlags fl )
    : QWidget( parent, fl )
{
	buildUI();
}

ProRataSearchPane::~ProRataSearchPane()
{
}

void ProRataSearchPane::buildUI()
{


	findLayout = new QHBoxLayout;
	findLayout->setSpacing( 0 );
	findLayout->setMargin( 0 );

	qleSearchString = new QLineEdit;
	findLayout->addWidget( qleSearchString );
	connect( qleSearchString, SIGNAL( returnPressed() ), this,   SLOT( searchClicked() ) );	

	qpbSearchButton = new QPushButton;
	connect( qpbSearchButton, SIGNAL( clicked() ), this, 
		SLOT( searchClicked() ) );
	qpbSearchButton->setText( tr( "Search" ) );
	findLayout->addWidget( qpbSearchButton );

	qpbShowAll = new QPushButton;
	connect( qpbShowAll, SIGNAL( clicked() ), this, 
		SLOT( showAllClicked() ) );
	qpbShowAll->setText( tr( "Show All" ) );
	findLayout->addWidget( qpbShowAll );

	setLayout( findLayout );
}

void ProRataSearchPane::searchClicked()
{
	QString sk( qleSearchString->text() );
	emit  searchStringSignal( qleSearchString->text() );
}

void ProRataSearchPane::showAllClicked()
{
	emit searchStringSignal( QString("") );
}

