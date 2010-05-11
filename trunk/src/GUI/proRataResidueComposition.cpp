
#include "proRataResidueComposition.h"

ProRataResidueComposition::ProRataResidueComposition( QWidget* qwParent, Qt::WFlags qwfFl )
	: QWidget( qwParent, qwfFl )
{
	priTreatment = NULL;
	priReference = NULL;
	buildUI();
	setValues();
}

ProRataResidueComposition::~ProRataResidueComposition()
{
}

void ProRataResidueComposition::buildUI()
{
	qgMainLayout = new QGridLayout;

	priTreatment = new ProRataIsotopologue;
	priTreatment->setName( "treatment" );

	priReference = new ProRataIsotopologue;
	priReference->setName( "reference" );

	qgMainLayout->addWidget( priTreatment, 0, 0 );
	qgMainLayout->addWidget( priReference, 1, 0 );

	setLayout( qgMainLayout );
}

void ProRataResidueComposition::setValues()
{

	QStringList qslTableData;
	qslTableData << "NTerm" <<  "CTerm" << "L" << "A" << "S" << "G" << "V" 
		   << "E" <<  "K" << "I" << "T" << "D" << "R" << "P" 
		   << "N" <<  "F" << "Q" << "Y" << "M" << "H" << "C" << "W";
	priTreatment->setData( 0, qslTableData );
	priReference->setData( 0, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "6" << "3" << "3" << "2" << "5" << "5" << "6" << "6" << "4" << "4" << "6" << "5" << "4" << "9" << "5" << "9" << "5" << "6" << "3" << "11";
	priTreatment->setData( 1, qslTableData );
	priReference->setData( 1, qslTableData );
	qslTableData.clear();

	qslTableData << "1"  << "1"  << "11" << "5"  << "5"  << "3"  << "9"  << "7"  << "12" << "11" << "7"  << "5"  << "12" << "7"  << "6"  << "9"  << "8"  << "9"  << "9"  << "7"  << "5"  << "10";
	priTreatment->setData( 2, qslTableData );
	priReference->setData( 2, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "1" << "1" << "1" << "2" << "1" << "1" << "3" << "1" << "1" << "2" << "3" << "1" << "1" << "2" << "1" << "2" << "2" << "1" << "1" << "1" << "1";
	priTreatment->setData( 3, qslTableData );
	priReference->setData( 3, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "1" << "1" << "1" << "1" << "1" << "1" << "2" << "1" << "1" << "1" << "4" << "1" << "2" << "1" << "2" << "1" << "1" << "3" << "1" << "2";
	priTreatment->setData( 4, qslTableData );
	priReference->setData( 4, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 5, qslTableData );
	priReference->setData( 5, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "1" << "0" << "1" << "0";
	priTreatment->setData( 6, qslTableData );
	priReference->setData( 6, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 7, qslTableData );
	priReference->setData( 7, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 8, qslTableData );
	priReference->setData( 8, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 9, qslTableData );
	priReference->setData( 9, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 10, qslTableData );
	priReference->setData( 10, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 11, qslTableData );
	priReference->setData( 11, qslTableData );
	qslTableData.clear();

	qslTableData << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0" << "0";
	priTreatment->setData( 12, qslTableData );
	priReference->setData( 12, qslTableData );
	qslTableData.clear();

}

