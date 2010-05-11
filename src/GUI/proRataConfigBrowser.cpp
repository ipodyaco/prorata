#include "proRataConfigBrowser.h"


ProRataConfigBrowser::ProRataConfigBrowser( QWidget* parent, Qt::WFlags fl )
	: QWidget( parent, fl ), qwParent( parent )
{
	bValidity = false;
	qsConfigFile = "";
	buildUI();
}

ProRataConfigBrowser::~ProRataConfigBrowser()
{

}

void ProRataConfigBrowser::buildUI()
{
	qvbMainLayout = new QVBoxLayout;

	qgbInputBox = new QGroupBox;
	qgbInputBox->setTitle( tr( "Configuration File" ) );

	qgbInputBoxLayout = new QGridLayout;
	qgbInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbInputBox->setLayout( qgbInputBoxLayout );

	qlFileLabel = new QLabel( qgbInputBox );
	qlFileLabel->setText( tr( "Name:" ) );

	qleFileEntry = new QLineEdit( qgbInputBox );
	qleFileEntry->setMinimumWidth( 100 );
	connect( qleFileEntry, SIGNAL( returnPressed() ),
			this, SLOT( validate() ) );
	connect( qleFileEntry, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validate() ) );

	qpbConfigFileBrowser = new QPushButton( qgbInputBox );
	qpbConfigFileBrowser->setText( tr( "&Browse..." ) );
	connect( qpbConfigFileBrowser, SIGNAL( clicked() ),
			this, SLOT( getConfigFile() ) );

	qpbNew = new QPushButton( qgbInputBox );
	qpbNew->setText( tr( "&New..." ) );
	connect( qpbNew, SIGNAL( clicked() ),
			this, SLOT( newConfigFile() ) );

	qgbInputBoxLayout->addWidget( qlFileLabel, 0, 0 );
	qgbInputBoxLayout->addWidget( qleFileEntry, 0, 1 );
	qgbInputBoxLayout->addWidget( qpbConfigFileBrowser, 0, 2 );
	qgbInputBoxLayout->addWidget( qpbNew, 0, 3 );

	qvbMainLayout->addWidget( qgbInputBox );

	setLayout( qvbMainLayout );
	setFixedHeight( sizeHint().height() );
}

void ProRataConfigBrowser::getConfigFile()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a ProRata configuration file", NULL,
			"Configuration files (*.xml)" );

	if ( !qsFName.isEmpty() )
	{       
		qleFileEntry->setText( qsFName );
		qsConfigFile = qsFName;
		bValidity = true;
	}       
}

void ProRataConfigBrowser::newConfigFile()
{
	ProRataConfigDialog *prcd = new ProRataConfigDialog( this, Qt::Dialog );
	prcd->setAttribute( Qt::WA_ShowModal );
	prcd->show();
}

bool ProRataConfigBrowser::isValid()
{
	return ( QFile::exists( qsConfigFile ) );
}

void ProRataConfigBrowser::validate()
{
	if ( !( QFile::exists( qleFileEntry->text() ) ) )
	{
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( false );
	}
	else
	{
		qsConfigFile = qleFileEntry->text();
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( true );
	}
}

