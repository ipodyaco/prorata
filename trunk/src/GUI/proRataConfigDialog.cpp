
#include "proRataConfigDialog.h"
#include <QMessageBox>

ProRataConfigDialog::ProRataConfigDialog( QWidget* qwParent, Qt::WFlags qwfFl )
	: QWidget( qwParent, qwfFl )
{
	setConfigVersion();
	buildUI();
	//setValues();
}

ProRataConfigDialog::~ProRataConfigDialog()
{

}

void ProRataConfigDialog::buildUI()
{
	setWindowTitle( tr("Configuration - ProRata") );
	qgMainLayout = new QGridLayout;

	prCfgSic = new ProRataCfgSIC;
	prQuant = new ProRataQuant;

	qtwMainTab = new QTabWidget;
	qtwMainTab->addTab( getIsotopologue( "reference" ), tr("Denominator Sample") );
	qtwMainTab->addTab( getIsotopologue( "treatment" ), tr("Numerator Sample") );
	qtwMainTab->addTab( prCfgSic, tr("Chromatogram") );
	qtwMainTab->addTab( prQuant, tr("Quantification") );

	qgMainLayout->addWidget( qtwMainTab, 0, 0 );

	qhbControlsLayout = new QHBoxLayout;
	qhbControlsLayout->setSpacing( 15 );

	qpbCancel = new QPushButton;
	qpbCancel->setText( tr( "&Cancel" ) );
	connect( qpbCancel, SIGNAL( clicked() ),
			this, SLOT( close() ) );
	qhbControlsLayout->addWidget( qpbCancel, Qt::AlignRight );

	qpbSaveAs = new QPushButton;
	qpbSaveAs->setText( tr( "Save..." ) );
	connect( qpbSaveAs, SIGNAL( clicked() ),
			this, SLOT( saveAs() ) );
	qhbControlsLayout->addWidget( qpbSaveAs, Qt::AlignRight );
/*
	qpbSave = new QPushButton;
	qpbSave->setText( tr( "&Save" ) );
	connect( qpbSave, SIGNAL( clicked() ),
			this, SLOT( save() ) );
	//qhbControlsLayout->addWidget( qpbSave, Qt::AlignRight );
*/
	qhbControlsLayout->insertStretch( 0, 10);

	qgMainLayout->addLayout( qhbControlsLayout, 1, 0 );

//	resize( 450, 350 );

	setLayout( qgMainLayout );

}

void ProRataConfigDialog::setValues()
{

}

QWidget*  ProRataConfigDialog::getIsotopologue( const QString & qsName )
{

	ProRataIsotopologue* priIsotop = new ProRataIsotopologue;
	priIsotop->setName( qsName );

	vector< int > viCurrentRow;
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );
	viCurrentRow.push_back( 0 );

//	<R>	NTerm,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 0;
	viCurrentRow[ 1 ] = 1;
	viCurrentRow[ 2 ] = 0;
	viCurrentRow[ 3 ] = 0;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "NTerm" ), viCurrentRow );


//	<R>	CTerm,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 0;
	viCurrentRow[ 1 ] = 1;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 0;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "CTerm" ), viCurrentRow );
	
//	<R>	L,	6,	11,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 6;
	viCurrentRow[ 1 ] = 11;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "L" ), viCurrentRow );
	
//	<R>	A,	3,	5,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 3;
	viCurrentRow[ 1 ] = 5;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "A" ), viCurrentRow );
	
//	<R>	S,	3,	5,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 3;
	viCurrentRow[ 1 ] = 5;
	viCurrentRow[ 2 ] = 2;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "S" ), viCurrentRow );
	
//	<R>	G,	2,	3,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 2;
	viCurrentRow[ 1 ] = 3;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "G" ), viCurrentRow );
	
//	<R>	V,	5,	9,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 5;
	viCurrentRow[ 1 ] = 9;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "V" ), viCurrentRow );
	
//	<R>	E,	5,	7,	3,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 5;
	viCurrentRow[ 1 ] = 7;
	viCurrentRow[ 2 ] = 3;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "E" ), viCurrentRow );
	
//	<R>	K,	6,	12,	1,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 6;
	viCurrentRow[ 1 ] = 12;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 2;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "K" ), viCurrentRow );
	
//	<R>	I,	6,	11,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 6;
	viCurrentRow[ 1 ] = 11;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "I" ), viCurrentRow );
	
//	<R>	T,	4,	7,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 4;
	viCurrentRow[ 1 ] = 7;
	viCurrentRow[ 2 ] = 2;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "T" ), viCurrentRow );
	
//	<R>	D,	4,	5,	3,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 4;
	viCurrentRow[ 1 ] = 5;
	viCurrentRow[ 2 ] = 3;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "D" ), viCurrentRow );
	
//	<R>	R,	6,	12,	1,	4,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 6;
	viCurrentRow[ 1 ] = 12;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 4;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "R" ), viCurrentRow );
	
//	<R>	P,	5,	7,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 5;
	viCurrentRow[ 1 ] = 7;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "P" ), viCurrentRow );
	
//	<R>	N,	4,	6,	2,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 4;
	viCurrentRow[ 1 ] = 6;
	viCurrentRow[ 2 ] = 2;
	viCurrentRow[ 3 ] = 2;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "N" ), viCurrentRow );
	
//	<R>	F,	9,	9,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 9;
	viCurrentRow[ 1 ] = 9;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "F" ), viCurrentRow );
	
//	<R>	Q,	5,	8,	2,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 5;
	viCurrentRow[ 1 ] = 8;
	viCurrentRow[ 2 ] = 2;
	viCurrentRow[ 3 ] = 2;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "Q" ), viCurrentRow );
	
//	<R>	Y,	9,	9,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 9;
	viCurrentRow[ 1 ] = 9;
	viCurrentRow[ 2 ] = 2;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "Y" ), viCurrentRow );
	
//	<R>	M,	5,	9,	1,	1,	0,	1,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 5;
	viCurrentRow[ 1 ] = 9;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 1;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "M" ), viCurrentRow );
	
//	<R>	H,	6,	7,	1,	3,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 6;
	viCurrentRow[ 1 ] = 7;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 3;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "H" ), viCurrentRow );
	
//	<R>	C,	3,	5,	1,	1,	0,	1,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 3;
	viCurrentRow[ 1 ] = 5;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 1;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 1;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "C" ), viCurrentRow );
	
//	<R>	W,	11,	10,	1,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>
	viCurrentRow[ 0 ] = 11;
	viCurrentRow[ 1 ] = 10;
	viCurrentRow[ 2 ] = 1;
	viCurrentRow[ 3 ] = 2;
	viCurrentRow[ 4 ] = 0;
	viCurrentRow[ 5 ] = 0;
	viCurrentRow[ 6 ] = 0;
	viCurrentRow[ 7 ] = 0;
	viCurrentRow[ 8 ] = 0;
	viCurrentRow[ 9 ] = 0;
	viCurrentRow[ 11 ] = 0;
	viCurrentRow[ 12 ] = 0;
	priIsotop->insertRow( QString( "W" ), viCurrentRow );
	
	return priIsotop;
}
/*
void ProRataConfigDialog::save()
{

#ifdef _WIN32
	QString qsFName = QDir::currentPath() + QString( "\\ProRataConfig.xml" );
#else
	QString qsFName = QDir::currentPath() + QString( "/ProRataConfig.xml" );
#endif


	if (QFile::exists(qsFName) &&
			QMessageBox::question(
				this,
				tr("Overwrite File? -- ProRata"),
				tr("A file called %1 already exists."
					"Do you want to overwrite it?")
				.arg(qsFName),
				tr("&Yes"), tr("&No"),
				QString(), 0, 1))
	{
		return;
	}

	buildXML( qsFName );

}
*/

void ProRataConfigDialog::saveAs()
{
/*
	QString qsFName = QFileDialog::getSaveFileName( this,
			tr( "Choose a file to save the settings."), 
			QDir::currentPath(), "XML Files (*.xml)", new QString( "ProRataConfig.xml" ) );
*/
	QString qsFName = QFileDialog::getSaveFileName( this,
			tr( "Choose a file to save the settings."), 
			QDir::homePath() + QDir::separator () + "ProRataConfig.xml", "XML Files (*.xml)" );
	if ( qsFName.isEmpty() )
	{       
		return;
	}       

	buildXML( qsFName );
	emit savedFilename( qsFName );
	close();

}

void ProRataConfigDialog::buildXML( const QString &qsFName )
{
	const int iIndentSize = 4;

	QFile qfNewConfigFile(qsFName);
	if (!qfNewConfigFile.open(QFile::WriteOnly | QFile::Text)) {
		QMessageBox::warning(this, tr("ProRata - Configuration"),
				tr("Cannot write file %1:\n%2.")
				.arg(qsFName)
				.arg(qfNewConfigFile.errorString()));
		return;
	}

	QTextStream qtsOut( &qfNewConfigFile );


	
	QDomProcessingInstruction qdpiProcess = qddXmlDoc.createProcessingInstruction(
			"xml", "version=\"1.0\"");
	qddXmlDoc.appendChild(qdpiProcess);
	

	//qddXmlDoc.clear();

	QDomElement qdeRoot = qddXmlDoc.createElement("CONFIG");
	qdeRoot.setAttribute( QString( "version" ), qsVersion );
	qddXmlDoc.appendChild( qdeRoot );

	///////////  SIC_EXTRACTION Section //////////////////
	//

	QStringList qslSicValues = prCfgSic->getValues();
	int iCtr = 0;

	QDomElement qdeSic = qddXmlDoc.createElement("SIC_EXTRACTION");
	qdeRoot.appendChild( qdeSic );
	addElement( qdeSic, "MS_FILE_TYPE", qslSicValues.at(iCtr++) );
	addElement( qdeSic, "ID_FILE_TYPE", qslSicValues.at(iCtr++) );


	QString qsPercentage = qslSicValues.at(iCtr++);

	//addElement( qdeSic, "FASTA_FILE", "................" );

	QDomElement qdeRetTimeInt = qddXmlDoc.createElement("RETENTION_TIME_INTERVAL");
	qdeSic.appendChild( qdeRetTimeInt );
	addElement( qdeRetTimeInt, "MINUTES_BEFORE_MS2", qslSicValues.at(iCtr++) );
	addElement( qdeRetTimeInt, "MINUTES_AFTER_MS2", qslSicValues.at(iCtr++) );
	addElement( qdeRetTimeInt, "MINUTES_BETWEEN_DUPLICATE_MS2", qslSicValues.at(iCtr++) );

	QDomElement qdeMass2CrgInt = qddXmlDoc.createElement("MASS_TO_CHARGE_INTERVAL");
	qdeSic.appendChild( qdeMass2CrgInt );
	addElement( qdeMass2CrgInt, "PLUS_MZ_ERROR", qslSicValues.at(iCtr++) );
	addElement( qdeMass2CrgInt, "MINUS_MZ_ERROR", qslSicValues.at(iCtr++) );
	addElement( qdeMass2CrgInt, "ISOTOPIC_ENVELOP_CUTOFF", qslSicValues.at(iCtr++) );

	QDomElement qdeAtomComp = qddXmlDoc.createElement("ATOM_ISOTOPIC_COMPOSITION");
	qdeSic.appendChild( qdeAtomComp );


	double dPercentage = qsPercentage.toDouble();
	dPercentage /= 100.0;

	qsPercentage = QString::number( dPercentage );

	QString qsRemPercentage( QString::number( 1 - dPercentage ) );

	QDomElement qdeC = qddXmlDoc.createElement("C");
	qdeAtomComp.appendChild( qdeC );
	addElement( qdeC, "MASS_DA", "	12.000000,	13.003355	" );
	addElement( qdeC, "NATURAL", "	0.9893,		0.0107		" );
	addElement( qdeC, "ENRICHED", "\t" + qsRemPercentage + ",\t\t" + qsPercentage + "\t\t" );

	QDomElement qdeH = qddXmlDoc.createElement("H");
	qdeAtomComp.appendChild( qdeH );
	addElement( qdeH, "MASS_DA", "	1.007825,	2.014102	" );
	addElement( qdeH, "NATURAL", "	0.999885,	0.000115	" );
	addElement( qdeH, "ENRICHED", "\t" +qsRemPercentage + ",\t\t" + qsPercentage + "\t\t" );

	QDomElement qdeO = qddXmlDoc.createElement("O");
	qdeAtomComp.appendChild( qdeO );
	addElement( qdeO, "MASS_DA", "	15.994915,	16.999132,	17.999160	" );
	addElement( qdeO, "NATURAL", "	0.99757,	0.00038,	0.00205		" );
	addElement( qdeO, "ENRICHED", "\t" + qsRemPercentage + ",\t\t0.0,\t\t" + qsPercentage  + "\t\t");

	QDomElement qdeN = qddXmlDoc.createElement("N");
	qdeAtomComp.appendChild( qdeN );
	addElement( qdeN, "MASS_DA", "	14.003074,	15.000109	" );
	addElement( qdeN, "NATURAL", "	0.99632,	0.00368		" );
	addElement( qdeN, "ENRICHED", "\t" + qsRemPercentage + ",\t\t" + qsPercentage  + "\t\t");

	QDomElement qdeP = qddXmlDoc.createElement("P");
	qdeAtomComp.appendChild( qdeP );
	addElement( qdeP, "MASS_DA", "	30.973762	" );
	addElement( qdeP, "NATURAL", "	1.0\t\t" );
	addElement( qdeP, "ENRICHED", "	1.0\t\t" );

	QDomElement qdeS = qddXmlDoc.createElement("S");
	qdeAtomComp.appendChild( qdeS );
	addElement( qdeS, "MASS_DA", "	31.972071,	32.971459,	33.967867,	35.967081	" );
	addElement( qdeS, "NATURAL", "	0.9493,		0.0076,		0.0429,		0.0002		" );
	addElement( qdeS, "ENRICHED", "\t" + qsRemPercentage + ",\t\t0.0,\t\t" + qsPercentage + ",\t\t0.0\t\t");

	QDomElement qdeResidue = qddXmlDoc.createElement("RESIDUE_ATOMIC_COMPOSITION");
	qdeSic.appendChild( qdeResidue );

	ProRataIsotopologue *priNum;
	ProRataIsotopologue *priDen;

	priNum = reinterpret_cast<ProRataIsotopologue *>(qtwMainTab->widget(1));

//	QStringList qslNumValues = priNum->getValues();
//	iCtr = 0;
//	int iRowSize = (qslNumValues.size()) - 1;

	int iColumnCount = 12;
	int i;
	int k;

	QStringList qsRowNameList;
	vector< vector< int > > vviTableContent;
	priNum->getData( qsRowNameList, vviTableContent );
	
	QDomElement qdeIsotop = qddXmlDoc.createElement("ISOTOPOLOGUE");
	qdeIsotop.setAttribute( QString( "name" ), priNum->getName() );
	qdeResidue.appendChild( qdeIsotop );

//	addComment( qdeIsotop, QString("	Name	C	H	O	N	P	S	C*	H*	O*	N*	P*	S*	") );
	for( i = 0; i < qsRowNameList.size(); i++ )
	{
		QString qsRowElements = QString("\t") + qsRowNameList[i];
		for( k = 0; k < iColumnCount; k++ )
		{
			qsRowElements = qsRowElements + QString(",\t") + QString::number( vviTableContent[i][k] );

		}
		addElement( qdeIsotop, "R", qsRowElements + "\t" );
	}

	priDen = reinterpret_cast<ProRataIsotopologue *>(qtwMainTab->widget(0));

	priDen->getData( qsRowNameList, vviTableContent );
	
//	QStringList qslDenValues = priDen->getValues();
//	iCtr = 0;
//	iRowSize = (qslDenValues.size()) - 1;

	QDomElement qdeIsotop2 = qddXmlDoc.createElement("ISOTOPOLOGUE");
	qdeIsotop2.setAttribute( QString( "name" ), priDen->getName() );
	qdeResidue.appendChild( qdeIsotop2 );

//	addComment( qdeIsotop2, "	Name	C	H	O	N	P	S	C*	H*	O*	N*	P*	S*	" );
	for( i = 0; i < qsRowNameList.size(); i++ )
	{
		QString qsRowElements = QString("\t") + qsRowNameList[i];
		for( k = 0; k < iColumnCount; k++ )
		{
			qsRowElements = qsRowElements + QString(",\t") + QString::number( vviTableContent[i][k] );

		}
		addElement( qdeIsotop2, "R", qsRowElements + "\t" );
	}

	////////////////////////////////////////////////////////////////



	///////////  PEPTIDE_QUANTIFICATION Section //////////////////
	
	QStringList qslQuantValues = prQuant->getValues();
	iCtr = 0;

	QDomElement qdePeptideQuant = qddXmlDoc.createElement("PEPTIDE_QUANTIFICATION");
	qdeRoot.appendChild( qdePeptideQuant );

	QDomElement qdePeakDect = qddXmlDoc.createElement("PEAK_DETECTION");
	qdePeptideQuant.appendChild( qdePeakDect );

	QDomElement qdeChroSmting = qddXmlDoc.createElement("CHROMATOGRAM_SMOOTHING");
	qdePeakDect.appendChild( qdeChroSmting );
	addElement( qdeChroSmting, "ORDER", qslQuantValues.at(iCtr++) );
	addElement( qdeChroSmting, "WINDOW_SIZE", qslQuantValues.at(iCtr++) );

	QDomElement qdePeakShift = qddXmlDoc.createElement("PEAK_SHIFT");
	qdePeakDect.appendChild( qdePeakShift );
	addElement( qdePeakShift, "LEFT", qslQuantValues.at(iCtr++) );
	addElement( qdePeakShift, "RIGHT", qslQuantValues.at(iCtr++) );

	QDomElement qdeAbunRt = qddXmlDoc.createElement("ABUNDANCE_RATIO");
	qdePeptideQuant.appendChild( qdeAbunRt );
	addElement( qdeAbunRt, "NUMERATOR_ISOTOPOLOGUE", priNum->getName() );
	addElement( qdeAbunRt, "DENOMINATOR_ISOTOPOLOGUE", priDen->getName() );
	
	QDomElement qdePepLog2Rt = qddXmlDoc.createElement("LOG2_RATIO");
	qdePeptideQuant.appendChild( qdePepLog2Rt );
	addElement( qdePepLog2Rt, "MINIMUM", qslQuantValues.at(iCtr++) );
	addElement( qdePepLog2Rt, "MAXIMUM", qslQuantValues.at(iCtr++) );

	addElement( qdePeptideQuant, "LOG2_SNR_CUTOFF", qslQuantValues.at(iCtr++) );

	addElement( qdePeptideQuant, "REMOVE_AMBIGUOUS_PEPTIDES", qslQuantValues.at(iCtr++) );
	////////////////////////////////////////////////////////////////

	///////////  PROTEIN_QUANTIFICATION Section //////////////////
	QDomElement qdeProteinQuant = qddXmlDoc.createElement("PROTEIN_QUANTIFICATION");
	qdeRoot.appendChild( qdeProteinQuant );
	addElement( qdeProteinQuant, "MIN_PEPTIDE_NUMBER", qslQuantValues.at(iCtr++) );
	addElement( qdeProteinQuant, "MAX_CI_WIDTH", qslQuantValues.at(iCtr++) );
	addElement( qdeProteinQuant, "MAX_LOG2_SNR", qslQuantValues.at(iCtr++) );

	QDomElement qdeLog2Ratio = qddXmlDoc.createElement("LOG2_RATIO");
	qdeProteinQuant.appendChild( qdeLog2Ratio );
	addElement( qdeLog2Ratio, "MINIMUM", qslQuantValues.at(iCtr++) );
	addElement( qdeLog2Ratio, "MAXIMUM", qslQuantValues.at(iCtr++) );

	addElement( qdeProteinQuant, "LOG2_RATIO_DISCRETIZATION", qslQuantValues.at(iCtr++) );

	addElement( qdeProteinQuant, "FASTA_FILE", QDir::convertSeparators(qslQuantValues.at(iCtr++)) );

	QDomElement qdeStdDev = qddXmlDoc.createElement("STANDARD_DEVIATION");
	qdeProteinQuant.appendChild( qdeStdDev );
	addElement( qdeStdDev, "SLOPE", "-0.26" );
	addElement( qdeStdDev, "INTERCEPT", "1.3" );

	QDomElement qdeMean = qddXmlDoc.createElement("MEAN");
	qdeProteinQuant.appendChild( qdeMean );
	addElement( qdeMean, "SLOPE", "1.2" );
	addElement( qdeMean, "INTERCEPT", "0" );

	addElement( qdeProteinQuant, "SMOOTHING_PROBABILITY_SPACE", "0.15" );
	addElement( qdeProteinQuant, "LN_LIKELIHOOD_CUTOFF_OFFSET", "1.96" );
	////////////////////////////////////////////////////////////////

	qddXmlDoc.save( qtsOut, iIndentSize);
	qddXmlDoc.clear();


}

void ProRataConfigDialog::addElement( QDomElement & qdeParent, const QString &qsTagName,
		const QString &qsValue )
{
	QDomElement qdeElement = qddXmlDoc.createElement( qsTagName );

	QDomText qdt = qddXmlDoc.createTextNode( qsValue );
	qdeElement.appendChild( qdt );

	qdeParent.appendChild( qdeElement );
}

void ProRataConfigDialog::addElement( QDomElement & qdeParent, const QString &qsTagName )
{
	QDomElement qdeElement = qddXmlDoc.createElement( qsTagName );
	qdeParent.appendChild( qdeElement );

}

void ProRataConfigDialog::addComment( QDomElement & qdeParent, const QString &qsValue )
{
	QDomComment qdcComment = qddXmlDoc.createComment( qsValue );
	qdeParent.appendChild( qdcComment );
//	qdeParent.insertBefore( qdcComment, NULL ); 

}
