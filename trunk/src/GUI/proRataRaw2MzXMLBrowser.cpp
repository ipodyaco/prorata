#include "proRataRaw2MzXMLBrowser.h"
#include <qglobal.h>


ProRataRaw2MzXMLBrowser::ProRataRaw2MzXMLBrowser( QWidget* parent )
	: QDialog( parent )
{
	bValidity = false;
	qsRawDataDirectory = "";
	qsOutputDataDirectory = "";
	buildUI();
}


ProRataRaw2MzXMLBrowser::ProRataRaw2MzXMLBrowser( QWidget* parent, const QString & qsDir )
	: QDialog( parent )
{
	bValidity = false;
	qsRawDataDirectory = "";
	qsOutputDataDirectory = qsDir;
	buildUI();
	if( qsDir != "" && qsDir != "/mzXML" )
		qgbOutputBox->setEnabled( false );
}

ProRataRaw2MzXMLBrowser::~ProRataRaw2MzXMLBrowser()
{

}

void ProRataRaw2MzXMLBrowser::buildUI()
{
	qvbClassLayout = new QVBoxLayout;

	QString qsDesc = "";

	qsDesc = "Select an input directory containing the raw files and\nan output directory to hold the mzXML files.Please make\ncertain that the instrument softwares have heen installed\nfor the given raw files.\n\nThe raw files are converted with the GNU General Public\nLicense programs ReAdW and Wolf. ReAdW and Wolf\nare not part of ProRata, but are co-distributed with\nProRata for your convenience.\n\nFor more information, please visit\nhttp://sashimi.sourceforge.net/software_glossolalia.html";

	//qlHeading = new QLabel;
	//qlHeading->setText( "Select an input directory containing the raw files and an output directory to hold the mzXML files.\nPlease make certain that the instrument softwares have heen installed for the given raw files.\nThe raw files are converted with the GNU General Public License programs ReAdW and Wolf.\nReAdW and Wolf are not part of ProRata, but are co-distributed with ProRata for your convenience.\nFor more information, please visit http://sashimi.sourceforge.net/software_glossolalia.html" );

	/*
	qfLine0 = new QFrame;
	qfLine0->setLineWidth( 1 );
	qfLine0->setFrameStyle( QFrame::HLine | QFrame::Sunken );
	*/

	this->setWhatsThis( qsDesc );

	
	qhbMainLayout = new QHBoxLayout;

	qgbConvProg = new QGroupBox;
	qgbConvProg->setTitle( tr( "Raw File Format" ) );

	qgbConvProgLayout = new QGridLayout;
	qgbConvProgLayout->setAlignment( Qt::AlignTop );

	qgbConvProg->setLayout( qgbConvProgLayout );

	qcbConvProg = new QComboBox;
	qcbConvProg->addItem("ThermoFinnigan Xcalibur");
	qcbConvProg->addItem("Micromass MassLynx");

	qgbConvProgLayout->addWidget( qcbConvProg, 0, 0 );

	qgbInputBox = new QGroupBox;
	qgbInputBox->setTitle( tr( "Raw File Directory" ) );

	qgbInputBoxLayout = new QGridLayout;
	qgbInputBoxLayout->setAlignment( Qt::AlignTop );

	qgbInputBox->setLayout( qgbInputBoxLayout );

	qlDirectoryLabel = new QLabel( qgbInputBox );
	qlDirectoryLabel->setText( tr( "Directory:" ) );

	qleDirectoryEntry = new QLineEdit( qgbInputBox );
	qleDirectoryEntry->setMinimumWidth( 250 );
	connect( qleDirectoryEntry, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validate() ) );

	qpbDirectoryBrowser = new QPushButton( qgbInputBox );
	qpbDirectoryBrowser->setText( tr( "&Browse..." ) );
	connect( qpbDirectoryBrowser, SIGNAL( clicked() ),
			this, SLOT( getRawDirectory() ) );

	qgbInputBoxLayout->addWidget( qlDirectoryLabel, 0, 0 );
	qgbInputBoxLayout->addWidget( qleDirectoryEntry, 0, 1 );
	qgbInputBoxLayout->addWidget( qpbDirectoryBrowser, 0, 2 );

	// Output
	//
	qgbOutputBox = new QGroupBox;
	qgbOutputBox->setTitle( tr( "Output mzXML Directory" ) );

	qgbOutputBoxLayout = new QGridLayout;
	qgbOutputBoxLayout->setAlignment( Qt::AlignTop );

	qgbOutputBox->setLayout( qgbOutputBoxLayout );

	qlOutputDirLabel = new QLabel( qgbOutputBox );
	qlOutputDirLabel->setText( tr( "Directory:" ) );

	qleOutputDirEntry = new QLineEdit( qgbOutputBox );
	qleOutputDirEntry->setMinimumWidth( 250 );
	qleOutputDirEntry->setText( qsOutputDataDirectory );
	connect( qleOutputDirEntry, SIGNAL( textChanged( const QString & ) ),
			this, SLOT( validateOuput() ) );

	qpbOutputDirectoryBrowser = new QPushButton( qgbOutputBox );
	qpbOutputDirectoryBrowser->setText( tr( "&Browse..." ) );
	connect( qpbOutputDirectoryBrowser, SIGNAL( clicked() ),
			this, SLOT( getOutputDirectory() ) );

	qgbOutputBoxLayout->addWidget( qlOutputDirLabel, 0, 0 );
	qgbOutputBoxLayout->addWidget( qleOutputDirEntry, 0, 1 );
	qgbOutputBoxLayout->addWidget( qpbOutputDirectoryBrowser, 0, 2 );
	// ////////////

	qhbMainLayout->addWidget( qgbConvProg );
	qhbMainLayout->addWidget( qgbInputBox );

	// Controls
	qfLine1 = new QFrame;
	qfLine1->setLineWidth( 1 );
	qfLine1->setFrameStyle( QFrame::HLine | QFrame::Sunken );

//	qpgbStatus = new QProgressBar;
//	qpgbStatus->setMinimumWidth( 300 );

	qpbCancelButton = new QPushButton(tr("Cancel"));
	connect( qpbCancelButton, SIGNAL( clicked() ),
			this, SLOT( close() ) );
	qpbOkButton = new QPushButton(tr("&Convert"));
	connect( qpbOkButton, SIGNAL( clicked() ),
			this, SLOT( convert() ) );

	qhblButtonLayout = new QHBoxLayout;
//	qhblButtonLayout->addWidget( qpgbStatus );
	qhblButtonLayout->addStretch(1);
	qhblButtonLayout->addWidget(qpbCancelButton);
	qhblButtonLayout->addWidget(qpbOkButton);

	/*
	qvbClassLayout->addWidget( qlHeading );
	qvbClassLayout->addWidget( qfLine0 );
	*/
	qvbClassLayout->addLayout( qhbMainLayout );
	qvbClassLayout->addWidget( qgbOutputBox );
	qvbClassLayout->addWidget( qfLine1 );
	qvbClassLayout->addLayout( qhblButtonLayout );

	setLayout( qvbClassLayout );
	setFixedHeight( sizeHint().height() );
	setWindowTitle( QString( "Raw2MzXML Conversion" ) );
}

void ProRataRaw2MzXMLBrowser::getRawDirectory()
{

	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose the directory that contains raw data files.",
			QDir::currentPath() );
	

	if ( !qsDirName.isEmpty() )
	{
		QDir::setCurrent( qsDirName );
		qleDirectoryEntry->setText( qsDirName );
		qsRawDataDirectory = qsDirName;
		bValidity = true;
	}
}

void ProRataRaw2MzXMLBrowser::getOutputDirectory()
{

	QString qsDirName = QFileDialog::getExistingDirectory(
			this,   
			"Choose the directory to store the output mzXML files.",
			QDir::currentPath() );

	if ( !qsDirName.isEmpty() )
	{
		QDir::setCurrent( qsDirName );
		qleOutputDirEntry->setText( qsDirName );
		qsOutputDataDirectory = qsDirName;
	}
}


bool ProRataRaw2MzXMLBrowser::isValid()
{
	QDir qdTemp( qsRawDataDirectory );
	return qdTemp.exists();
}

void ProRataRaw2MzXMLBrowser::validate()
{
	QDir qdTemp( qleDirectoryEntry->text() );

	if ( !( qdTemp.exists() ) )
	{
		qsRawDataDirectory = "";
		qsConvProgram = "";
	}
	else
	{
		qsRawDataDirectory = qleDirectoryEntry->text();
		qsConvProgram = qcbConvProg->currentText();
	}
}

void ProRataRaw2MzXMLBrowser::validateOuput()
{
	QDir qdTemp( qleOutputDirEntry->text() );

	if ( !( qdTemp.exists() ) )
	{
		qsOutputDataDirectory = "";
	}
	else
	{
		qsOutputDataDirectory = qleOutputDirEntry->text();
	}
}

void ProRataRaw2MzXMLBrowser::convert()
{
	// Output directory name.
	QDir qdOutput( qsOutputDataDirectory );
	if( !qdOutput.exists() )
	{
		if( !qdOutput.mkpath( qsOutputDataDirectory ) )
		{
			QMessageBox::information(this, "Error - ProRata", "Cannot make the output directory." );
		}
	}
	
	QDir qdRawDirectory( qsRawDataDirectory );
	QStringList qslFilters;
	qslFilters << "*.raw";
	QStringList qslFileList = qdRawDirectory.entryList( qslFilters, QDir::Files );

	
	int iMin = 1;
	int iMax = qslFileList.size();
//	qpgbStatus->setRange( iMin, iMax );
	QString qsExeFilename;


	QString qsCommand = "cmd.exe /c \"echo ";
	if( qsConvProgram == "ThermoFinnigan Xcalibur" )
	{
		qsCommand = qsCommand +  " ReAdW.exe centroid ";
	}
	else if( qsCommand == "Micromass MassLynx" )
	{
		qsCommand = qsCommand +  " wolf.exe ";
	}
	else
		qsCommand = " ";
	
	for (int i = iMin; i <= iMax; i++)
	{
		QString qsFilename = qslFileList.at(i-1);
		QString qsRawFile = " \"" + qsRawDataDirectory + "/" + qsFilename + "\" ";
		QString qsMZXMLfile = " \"" + qsOutputDataDirectory + "/" + qsFilename.replace( ".raw" , ".mzXML", Qt::CaseInsensitive) + "\" ";
		//QProcess qpProcess;
		if( qsConvProgram == "ThermoFinnigan Xcalibur" )
		{
			qsCommand = qsCommand + " && \"" + QCoreApplication::applicationDirPath() + "/ReAdW.exe\" "
				+ qsRawFile
				+ " c "
				+ qsMZXMLfile;
		}
		else if( qsCommand == "Micromass MassLynx" )
		{
			qsCommand =  qsCommand + " && \"" + QCoreApplication::applicationDirPath() + "/wolf.exe\" "
				+ qsRawFile
				+ " " + qsMZXMLfile;
		}
		else
			qsCommand = " ";

		
		//qDebug( "%s\n", qsCommand.toAscii().data() );
		//QMessageBox::information(this, "command ", qsCommand );	
		//qpProcess.execute( qsCommand );
		//QMessageBox::information(this, "output ", QString( qpProcess.readAllStandardOutput() ) );
		//QMessageBox::information(this, "error ", QString( qpProcess.readAllStandardError() ) ) ;

	//	qpgbStatus->setValue( i );
	}

	qsCommand = qsCommand + QString( " && pause \"");
	
//	system( qsCommand.toAscii().data() );
	QProcess::startDetached( qsCommand.toAscii().data() ); 
	emit savedDir( qsOutputDataDirectory );
	close();
}
