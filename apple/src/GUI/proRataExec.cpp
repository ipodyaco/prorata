#include "proRataExec.h"


ProRataExec::ProRataExec( QWidget* parent, Qt::WFlags fl )
	: QWidget( parent, fl )
{
	buildUI();
	setWindowFlags( Qt::Dialog );
}

ProRataExec::~ProRataExec()
{

}

void ProRataExec::buildUI()
{
	qvbMainLayout = new QVBoxLayout;

	qlHeading = new QLabel;
	qlHeading->setText( tr( "ProRata Execution Parameters" ) );

	qvbMainLayout->addWidget( qlHeading );

	qgbInputBox = new QGroupBox;
	qgbInputBox->setTitle( tr( "Input" ) );

	qgbInputBoxLayout = new QGridLayout;
	qgbInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbInputBox->setLayout( qgbInputBoxLayout );


	// Working directory stuff.
	qlWorkingDir = new QLabel( qgbInputBox );
	qlWorkingDir->setText( tr( "Working Directory" ) );

	qleWorkingDir = new QLineEdit( qgbInputBox );
	qleWorkingDir->setMinimumWidth( 200 );

	qpbWorkingDirBrowse = new QPushButton( qgbInputBox );
	qpbWorkingDirBrowse->setText( tr( "Browse..." ) );
	connect( qpbWorkingDirBrowse, SIGNAL( clicked() ),
			this, SLOT( workingDirBrowseSlot() ) );


	/*
	// Id file stuff.
	qlIdFile = new QLabel( qgbInputBox );
	qlIdFile->setText( tr( "Identification File" ) );

	qleIdFile = new QLineEdit( qgbInputBox );
	qleIdFile->setMinimumWidth( 200 );

	qpbIdFileBrowse = new QPushButton( qgbInputBox );
	qpbIdFileBrowse->setText( tr( "Browse..." ) );
	connect( qpbIdFileBrowse, SIGNAL( clicked() ),
			this, SLOT( idFileBrowseSlot() ) );
	*/

	// Add all the input components in a grid layout.
	qgbInputBoxLayout->addWidget( qlWorkingDir, 0, 0 );
	qgbInputBoxLayout->addWidget( qleWorkingDir, 0, 1 );
	qgbInputBoxLayout->addWidget( qpbWorkingDirBrowse, 0, 2 );

	/*
	qgbInputBoxLayout->addWidget( qlIdFile, 1, 0 );
	qgbInputBoxLayout->addWidget( qleIdFile, 1, 1 );
	qgbInputBoxLayout->addWidget( qpbIdFileBrowse, 1, 2 );
	*/

	qvbMainLayout->addWidget( qgbInputBox );

	qhblControls = new QHBoxLayout;
	spacer1 = new QSpacerItem( 51, 21, QSizePolicy::Expanding, QSizePolicy::Minimum );

	qpbCancel = new QPushButton;
	qpbCancel->setText( tr( "Cancel" ) );
	connect( qpbCancel, SIGNAL( clicked() ),
			this, SLOT( cancelSlot() ) );

	qpbExec = new QPushButton;
	qpbExec->setText( tr( "Execute" ) );
	connect( qpbExec, SIGNAL( clicked() ),
			this, SLOT( okSlot() ) );

	qhblControls->addItem( spacer1 );
	qhblControls->addWidget( qpbCancel );
	qhblControls->addWidget( qpbExec );

	qvbMainLayout->addLayout( qhblControls );

	setLayout( qvbMainLayout );

	setFixedHeight( sizeHint().height() );
	//resize( 450, 100 );
}

void ProRataExec::workingDirBrowseSlot()
{

	/*
	QString qsDirName = QFileDialog::getExistingDirectory(
			this,
			"Choose the working directory",
			".",
			QFileDialog::DontResolveSymlinks);
		*/
	
	QString qsDirName = QFileDialog::getExistingDirectory(
			this,
			"Choose the working directory",
			QDir::rootPath());

	if ( !qsDirName.isEmpty() )
	{       
		qleWorkingDir->setText( qsDirName );
	}       
}

void ProRataExec::idFileBrowseSlot()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a Identification file", NULL,
			"Id Files (*.txt)" );

	if ( !qsFName.isEmpty() )
	{       
		qleIdFile->setText( qsFName );
	}       
}

void ProRataExec::okSlot()
{
	QString qsDir = qleWorkingDir->text();
	if( qsDir == "" )
		return;
	QString qsDTASelect = qsDir + QString("//DTASelect-filter.txt");
	QString qsConfig = qsDir + QString( "//ProRataConfig.xml" );
	QString qsMZXML = qsDir + QString( "//mzXML" );
	QDir mzDir( qsMZXML );
	if( !QFile::exists( qsDTASelect ) )
	{
		QMessageBox::information( this, "Error", "DTASelect-filter.txt doesn't exists" );
		close();
		return;
	}
	if( !QFile::exists( qsConfig ) )
	{
		QMessageBox::information( this, "Error", "ProRataConfig.xml doesn't exists" );
		close();
		return;
	}
	if( !mzDir.exists() )
	{
		QMessageBox::information( this, "Error", "mzXML directory doesn't exists" );
		close();
		return;
	}
	
//	QProcess qpDirTest;
//	qpDirTest.start( "dir" );
	
	close();
}

void ProRataExec::cancelSlot()
{
	close();
}

