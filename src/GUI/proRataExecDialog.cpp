
#include "proRataExecDialog.h"
#include "proRataParameters.h"
#include "sicInfo.h"
#include "proteomeInfo.h"
//#include <fcntl.h>
//#include <iostream>

#include <QtAssistant/QAssistantClient>


using namespace std;

/*ProRataExecDialog::ProRataExecDialog( QWidget* qwParent, Qt::WFlags qwfFl )
	: QWidget( qwParent, qwfFl ) */
ProRataExecDialog::ProRataExecDialog( QWidget* qwParent )
	: QDialog( qwParent )
{
	qsWorkingDirectory = qsDTAResultsFile = qsCFGFile = 
	qsMZXMLDir = "";
	buildUI();
	setAttribute( Qt::WA_ShowModal );
	setModal( true );
}

ProRataExecDialog::~ProRataExecDialog()
{

}

void ProRataExecDialog::buildUI()
{
	qvbMainLayout = new QVBoxLayout;
/*
	qlHeading = new QLabel( tr("<center><b>ProRata Data Processing</b></center>" ) );
	qfLine1 = new QFrame;
	qfLine1->setLineWidth( 1 );
	qfLine1->setFrameStyle( QFrame::HLine | QFrame::Sunken );

	qlDescription = new QLabel( tr( "Please select a working directory to start data processing. The working directory should contain the input files: DTASelect-filter.txt, ProRataConfig.xml and a mzXML sub-directory storing all mzXML files. After the data processing is finished, open the result file ProRata_Quantification.qpr.xml in the working directory for data inspection." ) );
	qfLine2 = new QFrame;
	qfLine2->setLineWidth( 1 );
	qfLine2->setFrameStyle( QFrame::HLine | QFrame::Sunken );
*/

	QString qsDesc = "";


	// Working directory stuff.
	qsDesc = "Select a valid directory to store intermediate data (input and output)";

	qgbWDirInputBox = new QGroupBox;
	qgbWDirInputBox->setTitle( tr( "Working Directory" ) );
	//qgbWDirInputBox->setToolTip( QString( "Select a valid directory to store intermediate data (input and output)" ) );
	qgbWDirInputBox->setToolTip( qsDesc );
	qgbWDirInputBox->setWhatsThis( qsDesc );

	qgbWDirInputBoxLayout = new QGridLayout;
	qgbWDirInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbWDirInputBox->setLayout( qgbWDirInputBoxLayout );

	qlDirectoryLabel = new QLabel( qgbWDirInputBox );
	qlDirectoryLabel->setText( tr( "Directory:" ) );
	qlDirectoryLabel->setMinimumWidth( 60 );

	qleDirectoryEntry = new QLineEdit( qgbWDirInputBox );
	qleDirectoryEntry->setMinimumWidth( 300 );
	connect( qleDirectoryEntry, SIGNAL( textChanged( const QString & ) ), this, SLOT( validateWDir(const QString & ) ) );

	qpbDirectoryBrowser = new QPushButton( qgbWDirInputBox );
	qpbDirectoryBrowser->setText( tr( "&Browse..." ) );
	connect( qpbDirectoryBrowser, SIGNAL( clicked() ),
			this, SLOT( getWdirDirectory() ) );

	qgbWDirInputBoxLayout->addWidget( qlDirectoryLabel, 0, 0 );
	qgbWDirInputBoxLayout->addWidget( qleDirectoryEntry, 0, 1 );
	qgbWDirInputBoxLayout->addWidget( qpbDirectoryBrowser, 0, 2 );

	// DTA Select Stuff.
	qsDesc = "Select a valid Identification file (DTASelect-filter.txt)";

	qgbDTAInputBox = new QGroupBox;
	qgbDTAInputBox->setTitle( tr( "Identification File" ) );
	//qgbDTAInputBox->setToolTip( QString( "Select a valid Identification file (DTASelect-filter.txt)" ) );
	qgbDTAInputBox->setToolTip( qsDesc );
	qgbDTAInputBox->setWhatsThis( qsDesc );

	qgbDTAInputBoxLayout = new QGridLayout;
	qgbDTAInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbDTAInputBox->setLayout( qgbDTAInputBoxLayout );

	qlResultFile = new QLabel( qgbDTAInputBox );
	qlResultFile->setText( tr( "Filename:" ) );
	qlResultFile->setMinimumWidth( 60 );

	qleResultFile = new QLineEdit( qgbDTAInputBox );
	qleResultFile->setMinimumWidth( 300 );

	qpbResultFileBrowser = new QPushButton( qgbDTAInputBox );
	qpbResultFileBrowser->setText( tr( "&Browse..." ) );
	connect( qpbResultFileBrowser, SIGNAL( clicked() ), this, SLOT( getDTAResultsFile() ) );

	qpbCreatDTASelect = new QPushButton( qgbDTAInputBox );
	qpbCreatDTASelect->setText( tr( "&Create..." ) );
	connect( qpbCreatDTASelect, SIGNAL( clicked() ),
			this, SLOT( creatDTASelect() ) );
	qpbCreatDTASelect->setEnabled( false );

	qgbDTAInputBoxLayout->addWidget( qlResultFile, 0, 0 );
	qgbDTAInputBoxLayout->addWidget( qleResultFile, 0, 1 );
	qgbDTAInputBoxLayout->addWidget( qpbResultFileBrowser, 0, 2 );
	qgbDTAInputBoxLayout->addWidget( qpbCreatDTASelect, 0, 3 );
	
	// Configuration file Stuff.
	qsDesc = "Select a valid ProRata configuration file (ProRataConfig.xml)";

	qgbCFGInputBox = new QGroupBox;
	qgbCFGInputBox->setTitle( tr( "ProRata Configuration File" ) );
	qgbCFGInputBox->setToolTip( qsDesc );
	qgbCFGInputBox->setWhatsThis( qsDesc );

	qgbCFGInputBoxLayout = new QGridLayout;
	qgbCFGInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbCFGInputBox->setLayout( qgbCFGInputBoxLayout );

	qlFileLabel = new QLabel( qgbCFGInputBox );
	qlFileLabel->setText( tr( "Filename:" ) );
	qlFileLabel->setMinimumWidth( 60 );

	qleFileEntry = new QLineEdit( qgbCFGInputBox );
	//qleFileEntry->setMinimumWidth( 100 );
//	connect( qleFileEntry, SIGNAL( textChanged( const QString & ) ), this, SLOT( validateConfigFile(const QString &) ) );

	qpbConfigFileBrowser = new QPushButton( qgbCFGInputBox );
	qpbConfigFileBrowser->setText( tr( "&Browse..." ) );
	connect( qpbConfigFileBrowser, SIGNAL( clicked() ), this, SLOT( getConfigFile() ) );

	qpbNew = new QPushButton( qgbCFGInputBox );
	qpbNew->setText( tr( "&Create..." ) );
	connect( qpbNew, SIGNAL( clicked() ),
			this, SLOT( creatConfigFile() ) );

	qgbCFGInputBoxLayout->addWidget( qlFileLabel, 0, 0 );
	qgbCFGInputBoxLayout->addWidget( qleFileEntry, 0, 1 );
	qgbCFGInputBoxLayout->addWidget( qpbConfigFileBrowser, 0, 2 );
	qgbCFGInputBoxLayout->addWidget( qpbNew, 0, 3 );
	
	// MzXML Stuff.
	qsDesc = "Select a valid directory that has mzXML files";

	qgbMZXMLInputBox = new QGroupBox;
	qgbMZXMLInputBox->setTitle( tr( "mzXML Directory" ) );
	qgbMZXMLInputBox->setToolTip( qsDesc );
	qgbMZXMLInputBox->setWhatsThis( qsDesc );

	qgbMZXMLInputBoxLayout = new QGridLayout;
	qgbMZXMLInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbMZXMLInputBox->setLayout( qgbMZXMLInputBoxLayout );

	qlMZXMLDirectoryLabel = new QLabel( qgbMZXMLInputBox );
	qlMZXMLDirectoryLabel->setText( tr( "Directory:" ) );
	qlMZXMLDirectoryLabel->setMinimumWidth( 60 );

	qleMZXMLDirectoryEntry = new QLineEdit( qgbMZXMLInputBox );
	qleMZXMLDirectoryEntry->setMinimumWidth( 300 );
//	connect( qleMZXMLDirectoryEntry, SIGNAL( textChanged( const QString & ) ), this, SLOT( validateMZMLDir(const QString &) ) );

	qpbMZXMLDirectoryBrowser = new QPushButton( qgbMZXMLInputBox );
	qpbMZXMLDirectoryBrowser->setText( tr( "&Browse..." ) );
	connect( qpbMZXMLDirectoryBrowser, SIGNAL( clicked() ),
			this, SLOT( getMZXMLDirectory() ) );

	qpbMzXMLButton = new QPushButton(tr("&Create..."));
	qpbMzXMLButton->setToolTip( QString( "Click here to convert Raw data to mzXML files" ) );
	connect( qpbMzXMLButton, SIGNAL( clicked() ),
			this, SLOT( creatMzxml() ) );

	qgbMZXMLInputBoxLayout->addWidget( qlMZXMLDirectoryLabel, 0, 0 );
	qgbMZXMLInputBoxLayout->addWidget( qleMZXMLDirectoryEntry, 0, 1 );
	qgbMZXMLInputBoxLayout->addWidget( qpbMZXMLDirectoryBrowser, 0, 2 );
	qgbMZXMLInputBoxLayout->addWidget(qpbMzXMLButton, 0, 3);

	// Line 3

	qfLine3 = new QFrame;
	qfLine3->setLineWidth( 1 );
	qfLine3->setFrameStyle( QFrame::HLine | QFrame::Sunken );

	// Control Panel
	//
	qpbHelpButton = new QPushButton(tr("&Help"));
	connect( qpbHelpButton, SIGNAL( clicked() ), this, SLOT( helpExec() ) );

	qlStatus = new QLabel();
	qlStatus->setText(" ");

	qpbCancelButton = new QPushButton(tr("Cancel"));
	connect( qpbCancelButton, SIGNAL( clicked() ),
			this, SLOT( reject() ) );

	qpbProceedButton = new QPushButton(tr("&Start"));
//	qpbProceedButton->setEnabled(false);
	connect( qpbProceedButton, SIGNAL( clicked() ),	this, SLOT( proceed() ) );

	qhblButtonLayout = new QHBoxLayout;
	//qhblButtonLayout->addWidget(qpbHelpButton);
	qhblButtonLayout->addWidget(qlStatus);
	qhblButtonLayout->addStretch(1);
	qhblButtonLayout->addWidget(qpbCancelButton);
	qhblButtonLayout->addWidget(qpbProceedButton);

	

	//qvbMainLayout->addWidget( qlHeading );
//	qvbMainLayout->addWidget( qfLine1 );
//	qvbMainLayout->addWidget( qlDescription );
//	qvbMainLayout->addWidget( qfLine2 );
	qvbMainLayout->addWidget( qgbWDirInputBox );
	qvbMainLayout->addWidget( qgbDTAInputBox );
	qvbMainLayout->addWidget( qgbCFGInputBox );
	qvbMainLayout->addWidget( qgbMZXMLInputBox );
	qvbMainLayout->addWidget( qfLine3 );
	qvbMainLayout->addLayout(qhblButtonLayout);

	//qgbDTAInputBox->setEnabled( false );
	//qgbCFGInputBox->setEnabled( false );
	//qpbMzXMLButton->setEnabled( false );
	//qgbMZXMLInputBox->setEnabled( false );

	setLayout( qvbMainLayout );
	setFixedHeight( sizeHint().height() );
	setWindowTitle( QString( "Data Processing" ) );

}

void ProRataExecDialog::helpExec()
{
/*
	QMessageBox::information(this, "Help - ProRata", 
	"Please select a working directory and provide appropriate input data to start data processing with ProRata.\nThe working directory should contain DTASelect-filter.txt, ProRataConfig.xml and a mzXML sub-directory storing all mzXML files.\nThose inputs can be copied from elsewhere to the working diretory by clicking \"Browse...\" or be created by clicking \"Create...\"\nAfter the data processing is finished, open the output file ProRata_Quantification.qpr.xml to browse the quantification results." );
*/

	QString path = "C:\\cvs\\ProRata\\src\\GUI\\";
	//QAssistantClient *qacHelp = new QAssistantClient(path, this );
	QAssistantClient *qacHelp = new QAssistantClient(path );
	qacHelp->openAssistant();
	qacHelp->showPage(QString( "ProRataDataProcessing.html") );

	qacHelp->showPage( path + QString( "ProRataDataProcessing.html") );
	QMessageBox::information( this, "Help- temp", path + QString( "ProRataDataProcessing.html") );
}

void ProRataExecDialog::creatMzxml()
{
	ProRataRaw2MzXMLBrowser *prr2m = new ProRataRaw2MzXMLBrowser( this, qsWorkingDirectory + "/" + QString( "mzXML" ) );
	prr2m->setAttribute( Qt::WA_ShowModal );
	prr2m->show();
	connect( prr2m, SIGNAL( savedDir( const QString & ) ), qleMZXMLDirectoryEntry, SLOT( setText( const QString & )  ) );
}

void ProRataExecDialog::creatConfigFile()
{
	ProRataConfigDialog *prcd = new ProRataConfigDialog( this, Qt::Dialog );
	prcd->setAttribute( Qt::WA_ShowModal );
	prcd->show();
	connect( prcd, SIGNAL( savedFilename( const QString & ) ), qleFileEntry, SLOT( setText( const QString & )  ) );
}

void ProRataExecDialog::creatDTASelect()
{
	ProRataMerge *prDTASelectMerge = new ProRataMerge( this );
	prDTASelectMerge->setAttribute( Qt::WA_ShowModal );
	prDTASelectMerge->show();
}


void ProRataExecDialog::getWdirDirectory()
{
	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose the working directory for ProRata data processing.", QDir::rootPath() );

	if ( !qsDirName.isEmpty() )
	{
		qleDirectoryEntry->setText( qsDirName );
	}
}

void ProRataExecDialog::getDTAResultsFile()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a DTASelect-filter.txt", NULL,
			"DTASelect-filter (*.txt)" );

	if ( !qsFName.isEmpty() )
	{       
		qleResultFile->setText( qsFName );
	}       
}

void ProRataExecDialog::getConfigFile()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a ProRata configuration file", NULL,
			"ProRataConfig (*.xml)" );

	if ( !qsFName.isEmpty() )
	{       
		qleFileEntry->setText( qsFName );
	}       
}

void ProRataExecDialog::getMZXMLDirectory()
{

	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose the directory containing mzXML files.",
			 QDir::currentPath() ); 

	if ( !qsDirName.isEmpty() )
	{
		qleMZXMLDirectoryEntry->setText( qsDirName );
	}
}

bool ProRataExecDialog::validateWDir( const QString & qsDir)
{
	
//	QDir qdTemp( qleDirectoryEntry->text() );
	QDir qdTemp( qsDir );

	if ( qdTemp.exists() && qsDir != "" )
	{
		qsWorkingDirectory = qsDir;
		QDir::setCurrent( qsWorkingDirectory );
		//See if DTASelect_filter.txt is available.
		if ( QFile::exists( qleDirectoryEntry->text() + "/DTASelect-filter.txt" ) )
		{
			qsDTAResultsFile = qleDirectoryEntry->text() + "/DTASelect-filter.txt";
		}
		else
		{
			qsDTAResultsFile = "";
		}
		qleResultFile->setText( qsDTAResultsFile );			
		

		if ( QFile::exists( qleDirectoryEntry->text() + "/ProRataConfig.xml" ) )
		{
			qsCFGFile = qleDirectoryEntry->text() + "/ProRataConfig.xml";
		}
		else
		{
			qsCFGFile = "";
		}
		qleFileEntry->setText( qsCFGFile );

		QDir qdTemp2( qleDirectoryEntry->text() + "/mzXML");

		if ( qdTemp2.exists() )
		{
			qsMZXMLDir = qleDirectoryEntry->text() + "/mzXML";
		}
		else
		{
			qsMZXMLDir = "";
		}
		qleMZXMLDirectoryEntry->setText( qsMZXMLDir );

	}
	else
	{
		qsWorkingDirectory = "";		
		qsDTAResultsFile = "";
		qsCFGFile = "";
		qsMZXMLDir = "";
		qleResultFile->setText( "" );
		qleFileEntry->setText( "" );
		qleMZXMLDirectoryEntry->setText( "" );
	//	qgbDTAInputBox->setEnabled( false );
	//	qgbCFGInputBox->setEnabled( false );
	//	qgbMZXMLInputBox->setEnabled( false );
	}
	return true;

}

bool ProRataExecDialog::validateDTAFile()
{
	QString qsDTALocal = qsWorkingDirectory + "/DTASelect-filter.txt";
	QString qsDTAinput = qleResultFile->text();

	if( QFile::exists( qsDTAinput ) )
	{
		if( qsDTAinput == qsDTALocal )
		{
			return true;
		}
		else
		{
			if( !QFile::exists( qsDTALocal ) )
			{
				return QFile::copy( qsDTAinput, qsDTALocal );
			//	return true;
				
			}
			else
			{
				int iAns = QMessageBox::question( this,
					QString("DTASelect-filter.txt - ProRata"),
					"A DTASelect-fitler.txt file already exists in the working directory.\nDo you want to overwrite it with the given file\nor proceed with the existing DTASelect-filter.txt file?",	
					QString("&Overwrite"), 
					QString("&Proceed"),
					QString("&Cancel"), 0, 2);

				if( iAns == 0)
				{
					QFile::remove( qsDTALocal );
					return QFile::copy( qsDTAinput, qsDTALocal );
				//	return true;
				}
				else if( iAns == 1)
				{
					return true;
				}
				else
				{	
					return false;
				}
				
			}
		}

	}
	else
	{
		QString qsErrorMessage = "Please select a valid DTASelect-filter.txt file!";
		QMessageBox::information(this, "Error - ProRata ", qsErrorMessage );
		return false;	
	}
}

bool ProRataExecDialog::validateConfigFile()
{

	QString qsConfigLocal = qsWorkingDirectory + "/ProRataConfig.xml";
	QString qsConfigInput = qleFileEntry->text();

	if( QFile::exists( qsConfigInput ) )
	{
		if( qsConfigInput == qsConfigLocal )
		{
			return true;
		}
		else
		{
			if( !QFile::exists( qsConfigLocal ) )
			{
				return QFile::copy( qsConfigInput, qsConfigLocal );
			//	return true;
				
			}
			else
			{
				int iAns = QMessageBox::question( this,
					QString("ProRataConfig.xml - ProRata"),
					"A ProRataConfig.xml file already exists in the working directory.\nDo you want to overwrite it with the given file\nor proceed with the existing ProRataConfig.xml file?",	
					QString("&Overwrite"), 
					QString("&Proceed"),
					QString("&Cancel"), 0, 2);

				if( iAns == 0)
				{
					QFile::remove( qsConfigLocal );
					return QFile::copy( qsConfigInput, qsConfigLocal );
				//	return true;
				}
				else if( iAns == 1)
				{
					return true;
				}
				else
				{	
					return false;
				}
				
			}
		}

	}
	else
	{
		QString qsErrorMessage = "Please select a valid ProRataConfig.xml file!";
		QMessageBox::information(this, "ERROR ", qsErrorMessage );
		return false;	
	}

}

bool ProRataExecDialog::validateMZMLDir()
{

	QString qsMZdirLocal = qsWorkingDirectory + "/mzXML";
	QString qsMZdirInput = qleMZXMLDirectoryEntry->text();

	QDir qdDirLocal( qsMZdirLocal );
	QDir qdDirInput( qsMZdirInput );
	QDir qdWDir( qsWorkingDirectory );
	//QMessageBox::information(this, "Error - ProRata", qsMZdirInput );
	if( qdDirInput.exists() && qsMZdirInput != "" )
	{
		accept();
		if( qsMZdirInput == qsMZdirLocal )
		{
			return true;
		}
		else
		{
			/*int iAns = QMessageBox::question( this,
				QString("mzXML directory - ProRata"),
				QString("All mzXML files in the given directory will be copied to the ")+qsMZdirLocal+ 
				QString(" directory\nand overwrite the existing mzXML files with the same filenames."),	
				QString("&Proceed"), 
				QString("&Cancel"),
				QString(""), 0, 1);
				*/
				
			int iAns = 0;


			if( iAns == 1)
			{
				return false;
			}
			
			if( !qdDirLocal.exists() )
			{
				qdWDir.mkpath( qsMZdirLocal );	
			}
			QStringList qslFilters;
			qslFilters << "*.mzXML";
			QStringList qslFileList = qdDirInput.entryList( qslFilters, QDir::Files );
			for( int i = 0; i < qslFileList.size(); ++i )
			{
				QString qsFilename = qslFileList.at(i);
				QString qsInputFile = qsMZdirInput + "/" + qsFilename;
				QString qsLocalFile = qsMZdirLocal + "/" + qsFilename;
				//QFile::copy( qsInputFile, qsLocalFile );
				QFile qfInputFile;
				qfInputFile.setFileName( qsInputFile );
				qfInputFile.copy( qsLocalFile );
			}
			return true;

		}

	}
	else
	{
		QString qsErrorMessage = "Please select a valid directory containing mzXML files.";
		QMessageBox::information(this, "Error - ProRata", qsErrorMessage );
		return false;	
	}
}

void ProRataExecDialog::proceed()
{
	if ( qsWorkingDirectory == "" )
	{
		QMessageBox::information(this, "Error - ProRata", "Please select a valid working directory." );
		return;	
	}
	//qlStatus->setText( "Validating DTASelect File.." );
	qlStatus->setText( "" );

	if( !validateDTAFile() )
		return;

	qlStatus->setText( "" );
	
	if( !validateConfigFile() )
		return;

	//qlStatus->setText( "Validating MzXML Directory.." );
	qlStatus->setText( "" );

	if( !validateMZMLDir() )
		return;
	accept();

	if( QFile::exists( qsWorkingDirectory + "/ProRata_Quantification.qpr.xml" ) )
	{
		int iAns = QMessageBox::question(
				this,
				QString("QPR result file - ProRata"),
				"The file ProRata_Quantification.qpr.xml already exists in the working directory.\nDo you want to proceed to generate a new one or open the existing one?",
				QString("Proceed"), 
				QString("Open"),
				QString("Cancel"), 0, 2);

		if( iAns == 0)
		{
			// proceed
		}
		else if( iAns == 1)
		{
			emit qprFilename( qsWorkingDirectory + "/ProRata_Quantification.qpr.xml" );
			return;
		}
		else
		{	
			show();
			return;
		}
	}
	

	QDir qdXicDirectory(  qsWorkingDirectory + "/xic" );	

	QString qsPratio = "\"" + QCoreApplication::applicationDirPath() + QString( "/pratio.exe\" -w \"" ) 
		+ QDir::convertSeparators(qsWorkingDirectory) + QString( "\"" );
	QString qsSicForma = "\"" + QCoreApplication::applicationDirPath() + QString( "/sicForma.exe\" -w \"" )
	       	+ QDir::convertSeparators(qsWorkingDirectory) + QString( "\"" );
	QString qsSicFormaPRatioCommand = QString( "cmd.exe /c \" echo sicForma -w " ) + qsWorkingDirectory + QString( " && ") + qsSicForma
		+ QString( " && echo PRatio -w " ) + qsWorkingDirectory + QString( " && " ) + qsPratio
		+ QString( " && pause\"");
	QString qsPRatioCommand = QString( "cmd.exe /c \" echo PRatio -w " ) + qsWorkingDirectory + QString( " && ") + qsPratio
		+ QString( " && pause\"");
	
	if( !qdXicDirectory.exists() )
	{
		QProcess::startDetached( qsSicFormaPRatioCommand.toAscii().data() );

	}
	else
	{
		QString qsTextLabel = QString( "There already exists an xic directory in the working directory. \nDo you want to use the existing xic files or generate new xic files?") ;
		int iAns = QMessageBox::question(
				this,
				QString("ProRata -- sicForma."),
				qsTextLabel,
				QString("Use existing xic"), 
				QString("Generate new xic"),
				QString("Cancel"), 0, 2);

		if( iAns == 0)
		{
			QProcess::startDetached( qsPRatioCommand.toAscii().data() );

		}
		else if( iAns == 1)
		{

			QProcess::startDetached( qsSicFormaPRatioCommand.toAscii().data() );
		}
		else
		{	
			show();
			return;
		}
	}


	
	/*
	int fd = open("ProRataOutput.txt", O_WRONLY | O_CREAT | O_TRUNC, 0666);
	dup2( fd, 1 );
	dup2( fd, 2 );
	cout << " test dup2 1 " << endl;

	ProRataParameters params;
	if( !(params.setWorkingDirectory( qsWorkingDirectory.toAscii().data() ) ) )
		QMessageBox::information(this, "Problem ", "setWorkingDirectory	" );

	if( !(params.setIDFile( "DTASelect-filter.txt" ) ) )
		QMessageBox::information(this, "Problem ", "setIDFile	" );
	if( !(params.setConfigFile( "ProRataConfig.xml" ) ) )
		QMessageBox::information(this, "Problem ", "setConfigFile	" );


	SICinfo sicInfo;
	sicInfo.setFilename( params.getIDFilename() );
	sicInfo.process();
	ProteomeInfo mainProteomeInfo;
	if( !mainProteomeInfo.processPeptidesXIC() )
	{
		cout << "cannot process peptide XIC " << endl;
		QMessageBox::information(this, "Problem ", "mainProteomeInfo.processPeptidesXIC()" );
	}
	if( !mainProteomeInfo.processProteins() )
	{
		cout << "cannot process protein " << endl;
		QMessageBox::information(this, "Problem ", "mainProteomeInfo.processProteins()" );
	}
	mainProteomeInfo.writeFileQPR();
	mainProteomeInfo.writeFileTAB();
	*/

}

