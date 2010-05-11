#include "proRataDTASelect.h"


ProRataDTASelect::ProRataDTASelect( QWidget* parent, Qt::WFlags fl )
	: QWidget( parent, fl ), qwParent( parent )
{
	bValidity = false;
	qsResults = "";
	buildUI();
}

ProRataDTASelect::~ProRataDTASelect()
{

}

void ProRataDTASelect::buildUI()
{
	qvbMainLayout = new QVBoxLayout;

	qgbInputBox = new QGroupBox;
	qgbInputBox->setTitle( tr( "DTA Select Results" ) );

	qgbInputBoxLayout = new QGridLayout;
	qgbInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbInputBox->setLayout( qgbInputBoxLayout );

	qlResultFile = new QLabel( qgbInputBox );
	qlResultFile->setText( tr( "Results File:" ) );

	qleResultFile = new QLineEdit( qgbInputBox );
	qleResultFile->setMinimumWidth( 100 );
	connect( qleResultFile, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validate() ) );

	qpbResultFileBrowser = new QPushButton( qgbInputBox );
	qpbResultFileBrowser->setText( tr( "Browse..." ) );
	connect( qpbResultFileBrowser, SIGNAL( clicked() ),
			this, SLOT( getDTAResultsFile() ) );

	qgbInputBoxLayout->addWidget( qlResultFile, 0, 0 );
	qgbInputBoxLayout->addWidget( qleResultFile, 0, 1 );
	qgbInputBoxLayout->addWidget( qpbResultFileBrowser, 0, 2 );

	qvbMainLayout->addWidget( qgbInputBox );

	setLayout( qvbMainLayout );
	setFixedHeight( sizeHint().height() );
}

void ProRataDTASelect::getDTAResultsFile()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a result file", NULL,
			"DTA Select Results (*.txt)" );

	if ( !qsFName.isEmpty() )
	{       
		qleResultFile->setText( qsFName );
		qsResults = qsFName;
		bValidity = true;
	}       
}


void ProRataDTASelect::validate()
{
	QFile qfTemp( qleResultFile->text() );

	if ( !( qfTemp.exists() ) )
	{
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( false );
	}
	else
	{
		qsResults = qleResultFile->text();
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( true );
	}
}

