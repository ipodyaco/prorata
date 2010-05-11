#include "proRataWorkingDirBrowser.h"


ProRataWorkingDirBrowser::ProRataWorkingDirBrowser( QWidget* parent, Qt::WFlags fl )
	: QWidget( parent, fl ), qwParent( parent )
{
	bValidity = false;
	qsWorkingDirectory = "";
	buildUI();
}

ProRataWorkingDirBrowser::~ProRataWorkingDirBrowser()
{

}

void ProRataWorkingDirBrowser::buildUI()
{
	qhbMainLayout = new QHBoxLayout;

	qgbInputBox = new QGroupBox;
	qgbInputBox->setTitle( tr( "ProRata Working Directory" ) );

	qgbInputBoxLayout = new QGridLayout;
	qgbInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbInputBox->setLayout( qgbInputBoxLayout );

	qlDirectoryLabel = new QLabel( qgbInputBox );
	qlDirectoryLabel->setText( tr( "Name:" ) );

	qleDirectoryEntry = new QLineEdit( qgbInputBox );
	qleDirectoryEntry->setMinimumWidth( 100 );
	connect( qleDirectoryEntry, SIGNAL( returnPressed() ),
			this, SLOT( validate() ) );
	connect( qleDirectoryEntry, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validate() ) );

	qpbDirectoryBrowser = new QPushButton( qgbInputBox );
	qpbDirectoryBrowser->setText( tr( "&Browse..." ) );
	connect( qpbDirectoryBrowser, SIGNAL( clicked() ),
			this, SLOT( getDirectory() ) );

	qgbInputBoxLayout->addWidget( qlDirectoryLabel, 0, 0 );
	qgbInputBoxLayout->addWidget( qleDirectoryEntry, 0, 1 );
	qgbInputBoxLayout->addWidget( qpbDirectoryBrowser, 0, 2 );

	qhbMainLayout->addWidget( qgbInputBox );

	setLayout( qhbMainLayout );
	setFixedHeight( sizeHint().height() );
}

void ProRataWorkingDirBrowser::getDirectory()
{

	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose a working directory for ProRata pre-processing.",
			".",    
			QFileDialog::DontResolveSymlinks);

	if ( !qsDirName.isEmpty() )
	{
		qleDirectoryEntry->setText( qsDirName );
		qsWorkingDirectory = qsDirName;
		bValidity = true;
	}
}

bool ProRataWorkingDirBrowser::isValid()
{
	QDir qdTemp( qsWorkingDirectory );
	return qdTemp.exists();
}

void ProRataWorkingDirBrowser::validate()
{
	QDir qdTemp( qleDirectoryEntry->text() );

	if ( !( qdTemp.exists() ) )
	{
		/*
		QMessageBox::critical( this, "File not found! - ProRata",
				"File not found" );
				*/
		//qleDirectoryEntry->clear();
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( false );

	}
	else
	{
		qsWorkingDirectory = qleDirectoryEntry->text();
		(reinterpret_cast<ProRataPreProcessWizard*>(qwParent))->setButtonEnabled( true );
	}
}

