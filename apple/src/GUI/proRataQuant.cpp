
#include "proRataQuant.h"


ProRataQuant::ProRataQuant( QWidget* qwParent, Qt::WFlags qwfFl )
    : QWidget( qwParent, qwfFl )
{
	buildUI();
	setValues();
}

ProRataQuant::~ProRataQuant()
{

}

void ProRataQuant::setValues()
{


    qgbPeptide->setTitle( tr( "Peptide Quantification Parameters" ) );
    qtlSmOrd->setText( tr( "Smoothing Order:" ) );
    qtlLtPkSft->setText( tr( "Left Peak Shift:" ) );
    qtlMinLogRatio->setText( tr( "Minimum Log2 Ratio:" ) );
    qtlLogSNRCutff->setText( tr( "Log2 SNR Cutoff:" ) );
    
    qleSmOrd->clear();
    qleSmOrd->insertItem( 0, "2" );
    qleSmOrd->insertItem( 1, "3" );
    qleSmOrd->setCurrentIndex(0);
    
    qteLtPkSft->setRange(0, 50);
    qteLtPkSft->setSuffix( " Scans");
    qteLtPkSft->setValue( 0 );
    
    // minimum peptide log2ratio
    qleMinLogRatio->setRange(-10, -3);
    qleMinLogRatio->setValue( -7 );
    
    qleLogSNRCutff->setRange(0.0, 2.0);
    qleLogSNRCutff->setDecimals(1);
    qleLogSNRCutff->setSingleStep( 0.1 );
    qleLogSNRCutff->setValue(1.0);
    
    qtlSmWinSize->setText( tr( "Smoothing Window Size:" ) );
    qtlRtPkSft->setText( tr( "Right Peak Shift:" ) );
    qtlMaxLogRatio->setText( tr( "Maximum Log2 Ratio:" ) );
    qtlRmAmbPeps->setText( tr( "Remove Ambiguous Peptides:" ) );
    
    qleSmWinSize->clear();
    qleSmWinSize->insertItem( 0, "5" );
    qleSmWinSize->insertItem( 1, "7" );
    qleSmWinSize->insertItem( 2, "9" );
    qleSmWinSize->setCurrentIndex(1);
    
    qteRtPkSft->setRange(0, 50);
    qteRtPkSft->setSuffix( " scans");
    qteRtPkSft->setValue( 0 );
    
    qleMaxLogRatio->setRange(3, 10);
    qleMaxLogRatio->setValue( 7 );
    
    
    qcbRmAmbPeps->clear();
    qcbRmAmbPeps->insertItem( 0, tr( "True" ) );
    qcbRmAmbPeps->insertItem( 1, tr( "False" ) );
    
    qgbProtein->setTitle( tr( "Protein Quantification Parameters" ) );
//    qgbMinProtein->setTitle( tr( "" ) );
    qtlMinLogRatioPro->setText( tr( "Minimum Log2 Ratio:" ) );
    qtlMinPepNumber->setText( tr( "Minimum Peptide Number:" ) );
    
    qleMinLogRatioPro->setRange(-10, -3);
    qleMinLogRatioPro->setValue( -7 );
    
    qleMinPepNumber->setRange(1, 4);
    qleMinPepNumber->setSuffix( " peptides");
    qleMinPepNumber->setValue( 2 );
    
    qtlLogRatioDis->setText( tr( "Log2 Ratio Discretization:" ) );

    qleLogRatioDis->clear();
    qleLogRatioDis->insertItem( 0, "0.1" );
    qleLogRatioDis->insertItem( 1, "0.05" );
    qleLogRatioDis->insertItem( 2, "0.01" );
    qleLogRatioDis->setCurrentIndex(0);    

//    qgbProtein->setTitle( tr( "Protein" ) );
//    qgbMaxProtein->setTitle( tr( "" ) );
    qtlMaxLogRatioPro->setText( tr( "Maximum Log2 Ratio:" ) );
    qtlMaxLogSNR->setText( tr( "Maximum Log2 SNR:" ) );
    
    qleMaxLogRatioPro->setRange(3, 10);
    qleMaxLogRatioPro->setValue( 7 );
    
    qleMaxLogSNR->setRange(2.0, 6.0);
    qleMaxLogSNR->setDecimals(1);
    qleMaxLogSNR->setSingleStep( 0.1 );
    qleMaxLogSNR->setValue(3.0);

    qtlCIWidth->setText( tr( "Maximum CI Width:" ) );
    qleCIWidth->setRange(0, 15.0);
    qleCIWidth->setDecimals(1);
    qleCIWidth->setSingleStep( 0.5 );
    qleCIWidth->setValue(3.0);


}

void ProRataQuant::buildUI()
{

    qvbProRataQuantLayout = new QVBoxLayout;

    qgbPeptide = new QGroupBox;
    qgbPeptideLayout = new QHBoxLayout;
    qgbPeptideLayout->setAlignment( Qt::AlignTop );

    qgbPeptide->setLayout( qgbPeptideLayout );

    qvbPepLabels1 = new QVBoxLayout;

    qtlSmOrd = new QLabel( qgbPeptide );
    qvbPepLabels1->addWidget( qtlSmOrd );

    qtlLtPkSft = new QLabel( qgbPeptide );
    qvbPepLabels1->addWidget( qtlLtPkSft );

    qtlMinLogRatio = new QLabel( qgbPeptide );
    qvbPepLabels1->addWidget( qtlMinLogRatio );

    qtlLogSNRCutff = new QLabel( qgbPeptide );
    qvbPepLabels1->addWidget( qtlLogSNRCutff );

    qgbPeptideLayout->addLayout( qvbPepLabels1 );

    qvbPepEntry1 = new QVBoxLayout;

    qleSmOrd = new QComboBox( qgbPeptide );
    qvbPepEntry1->addWidget( qleSmOrd );

    qteLtPkSft = new QSpinBox( qgbPeptide );
    qvbPepEntry1->addWidget( qteLtPkSft );

    qleMinLogRatio = new QSpinBox( qgbPeptide );
    qvbPepEntry1->addWidget( qleMinLogRatio );

    qleLogSNRCutff = new QDoubleSpinBox( qgbPeptide );
    qvbPepEntry1->addWidget( qleLogSNRCutff );

    qgbPeptideLayout->addLayout( qvbPepEntry1 );

    qvbPepLabels2 = new QVBoxLayout;

    qtlSmWinSize = new QLabel( qgbPeptide );
    qvbPepLabels2->addWidget( qtlSmWinSize );

    qtlRtPkSft = new QLabel( qgbPeptide );
    qvbPepLabels2->addWidget( qtlRtPkSft );

    qtlMaxLogRatio = new QLabel( qgbPeptide );
    qvbPepLabels2->addWidget( qtlMaxLogRatio );

    qtlRmAmbPeps = new QLabel( qgbPeptide );
    qvbPepLabels2->addWidget( qtlRmAmbPeps );

    qgbPeptideLayout->addLayout( qvbPepLabels2 );

    qvbPepEntry2 = new QVBoxLayout;

    qleSmWinSize = new QComboBox( qgbPeptide );
    qvbPepEntry2->addWidget( qleSmWinSize );

    qteRtPkSft = new QSpinBox( qgbPeptide );
    qvbPepEntry2->addWidget( qteRtPkSft );

    qleMaxLogRatio = new QSpinBox( qgbPeptide );
    qvbPepEntry2->addWidget( qleMaxLogRatio );

    qcbRmAmbPeps = new QComboBox( qgbPeptide );
    qvbPepEntry2->addWidget( qcbRmAmbPeps );

    qgbPeptideLayout->addLayout( qvbPepEntry2 );

    qvbProRataQuantLayout->addWidget( qgbPeptide );
	    
    qgbProtein = new QGroupBox;
    //qgbProteinLayout = new QGridLayout( qgbProtein );
    //qgbProteinLayout->setAlignment( Qt::AlignTop );
    //qgbProtein->setLayout( qgbProteinLayout );

    // New additions //
    QGridLayout *qgbNewProteinLayout = new QGridLayout( qgbProtein );
    qgbProtein->setLayout( qgbNewProteinLayout );
    // Additions end here //

    /*

    qgbMinProtein = new QGroupBox( qgbProtein );
	qgbMinProtein->setFlat( true );
    qgbMinProteinLayout = new QHBoxLayout;
    qgbMinProteinLayout->setAlignment( Qt::AlignTop );

    qgbMinProtein->setLayout( qgbMinProteinLayout );
    */

    //qvbLogMinProLabels = new QVBoxLayout;

    qtlMinLogRatioPro = new QLabel;
    //qvbLogMinProLabels->addWidget( qtlMinLogRatioPro );

    qtlMinPepNumber = new QLabel;
    //qvbLogMinProLabels->addWidget( qtlMinPepNumber );

    //qgbMinProteinLayout->addLayout( qvbLogMinProLabels );

    //qvbLogMinProEntrys = new QVBoxLayout;
    
    qleMinLogRatioPro = new QSpinBox;
    //qvbLogMinProEntrys->addWidget( qleMinLogRatioPro );

    qleMinPepNumber = new QSpinBox;
    //qvbLogMinProEntrys->addWidget( qleMinPepNumber );

    //qgbMinProteinLayout->addLayout( qvbLogMinProEntrys );

    //qgbProteinLayout->addWidget( qgbMinProtein, 0, 0 );


    //qhbLgRtDis = new QHBoxLayout;
    qtlLogRatioDis = new QLabel( qgbProtein );
    //qhbLgRtDis->addWidget( qtlLogRatioDis );
    qleLogRatioDis = new QComboBox( qgbProtein );
    //qhbLgRtDis->addWidget( qleLogRatioDis );

    //qgbProteinLayout->addLayout( qhbLgRtDis, 1, 0 );

    /*
    qgbMaxProtein = new QGroupBox;
	qgbMaxProtein->setFlat( true );
    qgbMaxProteinLayout = new QVBoxLayout;
    qgbMaxProteinLayout->setAlignment( Qt::AlignTop );
    qgbMaxProtein->setLayout( qgbMaxProteinLayout );
    */
    
    //qvbLogMaxProLabels = new QVBoxLayout;
    qtlMaxLogRatioPro = new QLabel;
    //qvbLogMaxProLabels->addWidget( qtlMaxLogRatioPro );
    qtlMaxLogSNR = new QLabel;
    //qvbLogMaxProLabels->addWidget( qtlMaxLogSNR );


    //qvbLogMaxProEntrys = new QVBoxLayout;

    qleMaxLogRatioPro = new QSpinBox;
    //qvbLogMaxProEntrys->addWidget( qleMaxLogRatioPro );
    
    qleMaxLogSNR = new QDoubleSpinBox;
    //qvbLogMaxProEntrys->addWidget( qleMaxLogSNR );

    //qhbProMaxLogs = new QHBoxLayout;
    //qhbProMaxLogs->addLayout( qvbLogMaxProLabels );
    //qhbProMaxLogs->addLayout( qvbLogMaxProEntrys );

    //qgbMaxProteinLayout->addLayout( qhbProMaxLogs );

    //qhbCIWidth = new QHBoxLayout;
    qtlCIWidth = new QLabel;
    //qhbCIWidth->addWidget( qtlCIWidth );
    qleCIWidth = new QDoubleSpinBox;
    //qhbCIWidth->addWidget( qleCIWidth );

    //qgbMaxProteinLayout->addLayout( qhbCIWidth );

    //qgbProteinLayout->addWidget( qgbMaxProtein, 0, 1, 2, 1 );


    // New Additions //
    qgbNewProteinLayout->addWidget( qtlMinLogRatioPro, 0, 0 );
    qgbNewProteinLayout->addWidget( qleMinLogRatioPro, 0, 1 );
    qgbNewProteinLayout->addWidget( qtlMaxLogRatioPro, 0, 2 );
    qgbNewProteinLayout->addWidget( qleMaxLogRatioPro, 0, 3 );

    qgbNewProteinLayout->addWidget( qtlMinPepNumber, 1, 0 );
    qgbNewProteinLayout->addWidget( qleMinPepNumber, 1, 1 );
    qgbNewProteinLayout->addWidget( qtlMaxLogSNR, 1, 2 );
    qgbNewProteinLayout->addWidget( qleMaxLogSNR, 1, 3 );

    qgbNewProteinLayout->addWidget( qtlLogRatioDis, 2, 0 );
    qgbNewProteinLayout->addWidget( qleLogRatioDis, 2, 1 );
    qgbNewProteinLayout->addWidget( qtlCIWidth, 2, 2 );
    qgbNewProteinLayout->addWidget( qleCIWidth, 2, 3 );

    // End Additions //

    ///////////

    qgbFasta = new QGroupBox;
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

    //qgbProteinLayout->addWidget( qgbFasta, 2, 0, 1, 2 );

    ///////////////
    //
    // New Additions //
    qgbNewProteinLayout->addWidget( qgbFasta, 3, 0, 1, 4 );
    // End Additions //


    qvbProRataQuantLayout->addWidget( qgbProtein );
    qvbProRataQuantLayout->insertStretch(-1);

    setLayout( qvbProRataQuantLayout );

}


void ProRataQuant::fastaBrowseSlot()
{
	QString qsFName = QFileDialog::getOpenFileName( this,
			"Choose a FASTA file", NULL,
			"FASTA Files (*.fasta)" );

	if ( !qsFName.isEmpty() )
	{       
		qleFasta->setText( qsFName );
	}       
}

const QStringList &  ProRataQuant::getValues()
{
	qslValues.clear();

	// Dependency, Possible segfault if # of values are reduced.
	qslValues         /* Peptide values */
		<< qleSmOrd->currentText() << qleSmWinSize->currentText() << QString::number( qteLtPkSft->value() )
		<<  QString::number( qteRtPkSft->value() ) <<  QString::number( qleMinLogRatio->value() ) <<  QString::number( qleMaxLogRatio->value() ) 
		<< QString::number( qleLogSNRCutff->value() ) << qcbRmAmbPeps->currentText()
			/* Protein values */
		<< QString::number( qleMinPepNumber->value() ) << QString::number( qleCIWidth->value() ) << QString::number( qleMaxLogSNR->value() )
		<< QString::number( qleMinLogRatioPro->value() )  << QString::number( qleMaxLogRatioPro->value() )
		<< qleLogRatioDis->currentText() << qleFasta->text();

	return qslValues;
}

