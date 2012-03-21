#include "proRataCfgSICExt.h"


ProRataCfgSIC::ProRataCfgSIC( QWidget* parent, Qt::WFlags fl )
	: QWidget( parent, fl )
{
	buildUI();
	setValues();
}

ProRataCfgSIC::~ProRataCfgSIC()
{
}

void ProRataCfgSIC::buildUI()
{

	qvbMainLayout = new QVBoxLayout;


	qgbNewCfg = new QGroupBox;
	qgbNewCfg->setTitle( tr( "Selected Ion Chromatogram Extraction Configuration" ) );

	qgbNewCfgLayout = new QVBoxLayout;
	qgbNewCfgLayout->setAlignment( Qt::AlignTop );

	qgbNewCfg->setLayout( qgbNewCfgLayout );

	qlyMsIdFastaLayout = new QHBoxLayout;

	qgbMSFile = new QGroupBox( qgbNewCfg );
	qgbMSFile->setTitle( tr( "MS File Type" ) );

	qgbMSFileLayout = new QHBoxLayout;
	qgbMSFileLayout->setAlignment( Qt::AlignTop );

	qgbMSFile->setLayout( qgbMSFileLayout );

	qcbMSFile = new QComboBox( qgbMSFile );
	qcbMSFile->clear();
	qcbMSFile->insertItem( 0, tr( "mzXML" ) );
	qcbMSFile->setToolTip( "Select MS File type" );

	qgbMSFileLayout->addWidget( qcbMSFile );
	qlyMsIdFastaLayout->addWidget( qgbMSFile );

	qgbIDFile = new QGroupBox( qgbNewCfg );
	qgbIDFile->setTitle( tr( "ID File Type" ) );

	qgbIDFileLayout = new QHBoxLayout;
	qgbIDFileLayout->setAlignment( Qt::AlignTop );

	qgbIDFile->setLayout( qgbIDFileLayout );

	qcbIDFile = new QComboBox( qgbIDFile );
	qcbIDFile->clear();
	qcbIDFile->insertItem( 0, tr( "DTASelect" ) );
	qcbIDFile->setToolTip( "Select Identification file type" );

	qgbIDFileLayout->addWidget( qcbIDFile );
	qlyMsIdFastaLayout->addWidget( qgbIDFile );

	qgbConposition = new QGroupBox;
	qgbConposition->setTitle( tr( "Stable Isotope Enrichment" ) );

	qgbConpositionLayout = new QVBoxLayout;
	qgbConpositionLayout->setAlignment( Qt::AlignTop );

	qgbConposition->setLayout( qgbConpositionLayout );

	qlyEnrichmentLayout = new QHBoxLayout;
	qlEnrichmentPctg = new QLabel( qgbConposition );
	qlEnrichmentPctg->setText( tr("Enrichment Percentage") );
	qleEnrichmentPctg = new QDoubleSpinBox;
	qlyEnrichmentLayout->addWidget( qlEnrichmentPctg );
	qlyEnrichmentLayout->addWidget( qleEnrichmentPctg );

	qgbConpositionLayout->addLayout( qlyEnrichmentLayout );

	qlyMsIdFastaLayout->addWidget( qgbConposition );
	
	qgbNewCfgLayout->addLayout( qlyMsIdFastaLayout );

	

	/*
	qgbFasta = new QGroupBox( qgbNewCfg );
	qgbFasta->setTitle( tr( "FASTA File" ) );

	qgbFastaLayout = new QHBoxLayout;
	qgbFastaLayout->setAlignment( Qt::AlignTop );

	qgbFasta->setLayout( qgbFastaLayout );

	qleFasta = new QLineEdit( qgbFasta );
	qleFasta->setToolTip( "Give tha path of a FASTA file" );
	qgbFastaLayout->addWidget( qleFasta );

	qpbFasta = new QPushButton( qgbFasta );
	qpbFasta->setText( tr( "Bro&wse.." ) );

	connect( qpbFasta, SIGNAL( clicked() ), 
			this, SLOT( fastaBrowseSlot() ) );

	qpbFasta->setToolTip( "Select a FASTA file" );

	qgbFastaLayout->addWidget( qpbFasta );
	qlyMsIdFastaLayout->addWidget( qgbFasta );
	qgbNewCfgLayout->addLayout( qlyMsIdFastaLayout );
	*/

	qgbRetentionInt = new QGroupBox( qgbNewCfg );
	qgbRetentionInt->setTitle( tr( "Retention Time Interval [min]" ) );

	qgbRetentionIntLayout = new QHBoxLayout;
	qgbRetentionIntLayout->setAlignment( Qt::AlignTop );

	qgbRetentionInt->setLayout( qgbRetentionIntLayout );

	qvbRetentionLabelsLayout = new QVBoxLayout;
	qvbRetentionEntrysLayout = new QVBoxLayout;

	qlBeforeMS = new QLabel( qgbRetentionInt );
	qlBeforeMS->setText( tr( "Before MS2:" ) );
	//qgbRetentionIntLayout->addWidget( qlBeforeMS );
	qvbRetentionLabelsLayout->addWidget( qlBeforeMS );

	qleBeforeMS = new QDoubleSpinBox( qgbRetentionInt );
	qleBeforeMS->setToolTip( "Give retention time interval in minutes" );
	//qgbRetentionIntLayout->addWidget( qleBeforeMS );
	qvbRetentionEntrysLayout->addWidget( qleBeforeMS );

	qlAfterMS = new QLabel( qgbRetentionInt );
	qlAfterMS->setText( tr( "After MS2:" ) );
	//qgbRetentionIntLayout->addWidget( qlAfterMS );
	qvbRetentionLabelsLayout->addWidget( qlAfterMS );

	qleAfterMS = new QDoubleSpinBox( qgbRetentionInt );
	qleAfterMS->setToolTip( "Give retention time interval in minutes" );
	//qgbRetentionIntLayout->addWidget( qleAfterMS );
	qvbRetentionEntrysLayout->addWidget( qleAfterMS );

	qlBetMS = new QLabel( qgbRetentionInt );
	qlBetMS->setText( tr( "Between Duplicate MS2:" ) );
	//qgbRetentionIntLayout->addWidget( qlBetMS );
	qvbRetentionLabelsLayout->addWidget( qlBetMS );

	qleBetMS = new QDoubleSpinBox( qgbRetentionInt );
	qleBetMS->setToolTip( "Give retention time interval in minutes" );
	//qgbRetentionIntLayout->addWidget( qleBetMS );
	qvbRetentionEntrysLayout->addWidget( qleBetMS );

	qgbRetentionIntLayout->addLayout( qvbRetentionLabelsLayout );
	qgbRetentionIntLayout->addLayout( qvbRetentionEntrysLayout );
	//qgbRetentionIntLayout->insertSpacing( -1, 10 );
	//qgbRetentionIntLayout->addStretch();

	//qgbNewCfgLayout->addWidget( qgbRetentionInt );

	qgbMassToCharge = new QGroupBox;
	qgbMassToCharge->setTitle( tr( "Mass-to-Charge Isolation Window" ) );

	qgbMassToChargeLayout = new QHBoxLayout;
	qgbMassToChargeLayout->setAlignment( Qt::AlignTop );

	qgbMassToCharge->setLayout( qgbMassToChargeLayout );

	qvbMassToChargeLabelsLayout = new QVBoxLayout;
	qvbMassToChargeEntrysLayout = new QVBoxLayout;

	qlPlus = new QLabel( qgbMassToCharge );
	qlPlus->setText( tr( "Plus MZ Error:" ) );
	//qgbMassToChargeLayout->addWidget( qlPlus );
	qvbMassToChargeLabelsLayout->addWidget( qlPlus );

	qlePlus = new QDoubleSpinBox( qgbMassToCharge );
	qlePlus->setToolTip( "Give plus MZ error" );
	//qgbMassToChargeLayout->addWidget( qlePlus );
	qvbMassToChargeEntrysLayout->addWidget( qlePlus );

	qlMinus = new QLabel( qgbMassToCharge );
	qlMinus->setText( tr( "Minus MZ Error:" ) );
	//qgbMassToChargeLayout->addWidget( qlMinus );
	qvbMassToChargeLabelsLayout->addWidget( qlMinus );

	qleMinus = new QDoubleSpinBox( qgbMassToCharge );
//	qleMinus->setToolTip( "Give mass-to-charge interval" );
	//qgbMassToChargeLayout->addWidget( qleMinus );
	qvbMassToChargeEntrysLayout->addWidget( qleMinus );

	qlCutoff = new QLabel( qgbMassToCharge );
	qlCutoff->setText( tr( "Isotopic Abundance Cutoff:" ) );
	//qgbMassToChargeLayout->addWidget( qlCutoff );
	qvbMassToChargeLabelsLayout->addWidget( qlCutoff );

	qleCutff = new QDoubleSpinBox( qgbMassToCharge );
//	qleCutff->setToolTip( "Give mass-to-time interva" );
	//qgbMassToChargeLayout->addWidget( qleCutff );
	qvbMassToChargeEntrysLayout->addWidget( qleCutff );

	qgbMassToChargeLayout->addLayout( qvbMassToChargeLabelsLayout );
	qgbMassToChargeLayout->addLayout( qvbMassToChargeEntrysLayout );
	//qgbMassToChargeLayout->insertSpacing( -1, 10 );
	//qgbMassToChargeLayout->addStretch();

	// New Additions //
	
	QHBoxLayout * qhbNewLayout = new QHBoxLayout;
	qhbNewLayout->addWidget( qgbRetentionInt );
	qhbNewLayout->addWidget( qgbMassToCharge );
	qgbNewCfgLayout->addLayout( qhbNewLayout );

	// End Additions //
	//qgbNewCfgLayout->addWidget( qgbMassToCharge );


	qvbMainLayout->addWidget( qgbNewCfg );
	//qvbMainLayout->addSpacing( 10 );

	setLayout( qvbMainLayout );

}

/*
void ProRataCfgSIC::existingBrowseSlot()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose an existing file", NULL,
			"ProRataConfig Files (*.xml)" );

	if ( !qsFName.isEmpty() )
	{       
		qleExistingFile->setText( qsFName );
	}       
	qgbNewCfg->setChecked( FALSE );
}
*/

/*
void ProRataCfgSIC::fastaBrowseSlot()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a FASTA file", NULL,
			"FASTA Files (*.fasta)" );

	if ( !qsFName.isEmpty() )
	{       
		qleFasta->setText( qsFName );
	}       
}
*/

/*
void ProRataCfgSIC::addCompositionSlot()
{
	QString qsNewAtom = "";

	qsNewAtom = qteComposition->toPlainText();
	qsNewAtom += "\n";
	qsNewAtom = qsNewAtom + convertToXML( 
			convertToXML( qleMassDA->text(), QString( "MASS_DA" ) ) + 
			convertToXML( qleNatural->text(), QString( "NATURAL" ) ) + 
			convertToXML( qleEnriched->text(), QString( "ENRICHED" ) ), 
			qcbAtom->currentText() );
	qteComposition->setPlainText( qsNewAtom );
}

QString ProRataCfgSIC::convertToXML( QString qsContent,
		QString qsTag )
{
	QString qsNewStr;
	qsNewStr = "<" + qsTag + ">" + qsContent +
		"</" + qsTag + ">";

	return qsNewStr;
}
*/

void ProRataCfgSIC::setValues()
{
	qleBeforeMS->setRange(1.0, 10.0);
	qleBeforeMS->setDecimals(1);
	qleBeforeMS->setSuffix( " min");
	qleBeforeMS->setSingleStep( 0.5 );
	qleBeforeMS->setValue(2.0);	


	qleAfterMS->setRange(1.0, 10.0);
	qleAfterMS->setDecimals(1);
	qleAfterMS->setSuffix( " min");
	qleAfterMS->setSingleStep( 0.5 );
	qleAfterMS->setValue( 2.0);	

	qleBetMS->setRange(0.5, 10.0);
	qleBetMS->setDecimals(1);
	qleBetMS->setSuffix( " min");
	qleBetMS->setSingleStep( 0.5 );
	qleBetMS->setValue(2.0);	

	qlePlus->setRange(0.0, 2.0);
	qlePlus->setDecimals(2);
	qlePlus->setSuffix( " mu");
	qlePlus->setSingleStep( 0.1 );
	qlePlus->setValue(0.5);	
	
	qleMinus->setRange(0.0, 2.0);
	qleMinus->setDecimals(2);
	qleMinus->setSuffix( " mu");
	qleMinus->setSingleStep( 0.1 );
	qleMinus->setValue(0.5);	
	
	qleCutff->setRange(0.0, 100.0);
	qleCutff->setDecimals(1);
	qleCutff->setSuffix( " %");
	qleCutff->setSingleStep( 1 );
	qleCutff->setValue( 10 );	

	qleEnrichmentPctg->setRange(0.0, 100.0);
	qleEnrichmentPctg->setDecimals(1);
	qleEnrichmentPctg->setSuffix( " %");
	qleEnrichmentPctg->setSingleStep( 0.5 );
	qleEnrichmentPctg->setValue(98.0);	
}

const QStringList &  ProRataCfgSIC::getValues()
{
	qslValues.clear();

	qslValues <<   qcbMSFile->currentText() << qcbIDFile->currentText()
		<< QString::number( qleEnrichmentPctg->value() ) 
		<< QString::number( qleBeforeMS->value() ) 
	       	<<  QString::number( qleAfterMS->value() )  
		<<  QString::number( qleBetMS->value() )  
		<<  QString::number( qlePlus->value() ) 
		<<  QString::number( qleMinus->value() ) 
		<<  QString::number( (qleCutff->value()/100.0) ); 

	return qslValues;
}
