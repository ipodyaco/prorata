
#include "proRataMerge.h"

ProRataMerge::ProRataMerge( QWidget* qwParent )
	: QDialog( qwParent )
/*	: QWidget( qwParent, qwfFl ) */
{
	buildUI();
	setValues();

	qsOutputDirectory = "";
	setAttribute( Qt::WA_ShowModal );
	setModal( true );
}

ProRataMerge::~ProRataMerge()
{
}

void ProRataMerge::buildUI()
{
	qvbProRataMergeLayout = new QVBoxLayout;

	qbgInput = new QGroupBox;
	qbgInputLayout = new QGridLayout;

	qbgInput->setLayout( qbgInputLayout );

	qlwInput = new QListWidget;
	qbgInputLayout->addWidget( qlwInput, 0, 0, 2, 1 );

	qpbAddDir = new QPushButton;
	connect( qpbAddDir, SIGNAL( clicked() ),
			this, SLOT( addDirectory() ) );

	qbgInputLayout->addWidget( qpbAddDir, 0, 1 );

	qpbRemoveDir = new QPushButton;
	connect( qpbRemoveDir, SIGNAL( clicked() ),
			this, SLOT( removeDirectory() ) );
	qbgInputLayout->addWidget( qpbRemoveDir, 1, 1 );

	qvbProRataMergeLayout->addWidget( qbgInput );

	qbgOutput = new QGroupBox;
	qbgOutputLayout = new QHBoxLayout;

	qbgOutput->setLayout( qbgOutputLayout );

	qleOutput = new QLineEdit;
	connect( qleOutput, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validateOutputDir(const QString &) ) );

	qbgOutputLayout->addWidget( qleOutput );

	qpbOuputBrowse = new QPushButton;
	connect( qpbOuputBrowse, SIGNAL( clicked() ),
			this, SLOT( browseForOutputDirectory() ) );
	qbgOutputLayout->addWidget( qpbOuputBrowse );

	qvbProRataMergeLayout->addWidget( qbgOutput );

	// Separator
	qfLine1 = new QFrame;
	qfLine1->setLineWidth( 1 );
	qfLine1->setFrameStyle( QFrame::HLine | QFrame::Sunken );

	qvbProRataMergeLayout->addWidget( qfLine1 );

	// Controls

	qpbCancelButton = new QPushButton(tr("Cancel"));
	connect( qpbCancelButton, SIGNAL( clicked() ),
			this, SLOT( close() ) );
	qpbMergeButton = new QPushButton(tr("&Merge"));
	qpbMergeButton->setEnabled(false);

	connect( qpbMergeButton, SIGNAL( clicked() ),
			this, SLOT( merge() ) );

	qhblButtonLayout = new QHBoxLayout;
	qhblButtonLayout->addStretch(1);
	qhblButtonLayout->addWidget(qpbCancelButton);
	qhblButtonLayout->addWidget(qpbMergeButton);

	qvbProRataMergeLayout->addLayout( qhblButtonLayout );

	setLayout( qvbProRataMergeLayout );
}

void ProRataMerge::setValues()
{
	qbgInput->setTitle( tr( "Input Directories" ) );
	qlwInput->clear();
	qpbAddDir->setText( tr( "Add Directory.." ) );
	qpbRemoveDir->setText( tr( "Remove Directory" ) );
	qbgOutput->setTitle( tr( "Output Directory" ) );
	qpbOuputBrowse->setText( tr( "Browse.." ) );

	setWindowTitle( tr( "Merge Directories - ProRata" ) );
	setFixedHeight( sizeHint().height() );

}

void ProRataMerge::merge()
{
	accept();

	// Get the list of files.
	for (int i = 0; i < qlwInput->count() ; i++ )
	{
		QString qsDir = (qlwInput->item( i ))->text();
		QDir dir( qsDir );
		QString qsBaseName = qsDir.section( '/', -1 );
		dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);

		QFileInfoList list = dir.entryInfoList();
		for (int j = 0; j < list.size(); ++j) {
			QFileInfo fileInfo = list.at(j);
			QString qsNewFile = qsOutputDirectory + "/" + qsBaseName + "_" + fileInfo.fileName();
			//QFile::copy( fileInfo.absoluteFilePath(), qsNewFile );
			QFile qfInputFile;
			qfInputFile.setFileName( fileInfo.absoluteFilePath() );
			qfInputFile.copy( qsNewFile );
		}

		dir.setFilter(QDir::AllDirs | QDir::NoSymLinks);

		QFileInfoList qfileDirList = dir.entryInfoList();
		QString qsDirList = "";
		for (int j = 0; j < qfileDirList.size(); ++j) 
		{
			QFileInfo dirInfo = qfileDirList.at(j);

			if ( dirInfo.fileName() == "." || dirInfo.fileName() == ".." )
			{
				continue;
			}
			QDir qdCurrentPointer;
			QString qsDestDirectory = qsOutputDirectory + "/" + qsBaseName + "_" + dirInfo.fileName();
			qdCurrentPointer.mkdir( qsDestDirectory );

			QDir sourceDir( dirInfo.absoluteFilePath() );
			sourceDir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);

			QFileInfoList list = sourceDir.entryInfoList();
			for (int j = 0; j < list.size(); ++j) {
				QFileInfo dirInfo = list.at(j);
				QString qsNewFile = qsDestDirectory + "/" + qsBaseName + "_" + dirInfo.fileName();
				//QFile::copy( dirInfo.absoluteFilePath(), qsNewFile );
				QFile qfInputFile;
				qfInputFile.setFileName( dirInfo.absoluteFilePath() );
				qfInputFile.copy( qsNewFile );
			}

		}
		//QMessageBox::information(this, "dir list ", qsDirList );


	}
}

void ProRataMerge::addDirectory()
{
	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose a directory to merge.",
			QDir::rootPath());

	if ( !qsDirName.isEmpty() )
	{
		qlwInput->addItem( new QListWidgetItem( qsDirName ) );
	}

}

void ProRataMerge::removeDirectory()
{
	qlwInput->takeItem( qlwInput->currentRow() );
}

void ProRataMerge::browseForOutputDirectory()
{
	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose a directory to save the merge output.",
			QDir::rootPath());

	if ( !qsDirName.isEmpty() )
	{
		qleOutput->setText( qsDirName );
	}
}

void ProRataMerge::validateOutputDir( const QString & qsDir)
{

	QDir qdTemp2( qsDir );

	if ( qdTemp2.exists() )
	{
		qpbMergeButton->setEnabled(true);
		qsOutputDirectory = qleOutput->text();
	}
	else
	{
		qsOutputDirectory = "";
	}
}
