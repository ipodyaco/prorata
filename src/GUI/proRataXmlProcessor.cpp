
#include "proRataXmlProcessor.h"

ProRataXmlProcessors::ProRataXmlProcessors( const QString & qsFName )
{
	qsFilename = qsFName;

	prptProteinSeletedData = new ProRataProteinTuple();
	prdtPeptideSeletedData = new ProRataPeptideTuple();


	if ( ! processFile( qsFilename ) )
	{
		cout << "Problem processing xml file" << endl;
	}

	qdeCurrentProteinElement = qdeFirstProteinElement;
	qdeCurrentPeptideElement = qdeFirstPeptideElement;

}

ProRataXmlProcessors::ProRataXmlProcessors()
{
	prptProteinSeletedData = new ProRataProteinTuple();
	prdtPeptideSeletedData = new ProRataPeptideTuple();
}

void ProRataXmlProcessors::setFilename( const QString & qsFName )
{
	qsFilename = qsFName;

	if ( ! processFile( qsFilename ) )
	{
		cout << "Problem processing xml file" << endl;
	}

	qdeCurrentProteinElement = qdeFirstProteinElement;
	qdeCurrentPeptideElement = qdeFirstPeptideElement;

}

ProRataXmlProcessors::~ProRataXmlProcessors()
{
	delete prptProteinSeletedData;
	delete prdtPeptideSeletedData;
}

bool ProRataXmlProcessors::processProtein( QDomElement &qdeProtein )
{
	if ( ! qdeProtein.isNull() )
	{
		qdeFirstPeptideElement = qdeProtein.firstChildElement( "CHROMATOGRAM" );

		prptProteinSeletedData->setLocus( getText( qdeProtein, QString( "LOCUS" ) ) );
		prptProteinSeletedData->setDescription( getText( qdeProtein, QString( "DESCRIPTION" ) ) );
		prptProteinSeletedData->setProteinLog2Ratio( getText( qdeProtein, QString( "LOG2_RATIO" ) ).toDouble() );
		prptProteinSeletedData->setCIWidth( getText( qdeProtein, QString( "UPPER_CONFIDENCE_INTERVAL" ) ).toDouble() -
				getText( qdeProtein, QString( "UPPER_CONFIDENCE_INTERVAL" ) ).toDouble() );

		return true;
	}
	else
	{
		return false;
	}
}

bool ProRataXmlProcessors::processPeptide( QDomElement & qdePeptide )
{
	if ( ! qdePeptide.isNull() )
	{
		prdtPeptideSeletedData->setChroFilename( getText( qdePeptide, QString( "CHROMATOGRAM_FILENAME" ) ) );
		prdtPeptideSeletedData->setChroIdentifier( getText( qdePeptide, QString( "CHROMATOGRAM_IDENTIFIER" ) ) );
		prdtPeptideSeletedData->setChroValidity( getText( qdePeptide, QString( "VALIDITY" ) ) );
		prdtPeptideSeletedData->setChroSequence( getText( qdePeptide, QString( "SEQUENCE" ) ) );
		prdtPeptideSeletedData->setChroLog2Ratio( getText( qdePeptide, QString( "LOG2_RATIO" ) ).toDouble() );
		prdtPeptideSeletedData->setEigenvalueRatio( getText( qdePeptide, QString( "EIGENVALUE_RATIO" ) ).toDouble() );
		prdtPeptideSeletedData->setNTerminalDistance( getText( qdePeptide, QString( "N_TERMINAL_DISTANCE" ) ).toInt() );
		prdtPeptideSeletedData->setCTerminalDistance( getText( qdePeptide, QString( "C_TERMINAL_DISTANCE" ) ).toInt() );
		prdtPeptideSeletedData->setChargeState( getText( qdePeptide, QString( "CHARGE_STATE" ) ).toInt() );

		return true;
	}
	else
	{
		return false;
	}

}

ProRataProteinTuple * ProRataXmlProcessors::getProtein( const QString & qsLocus )
{
	setProtein( qsLocus );

	if ( ! qdeCurrentProteinElement.isNull() )
	{
		if ( processProtein( qdeCurrentProteinElement ) )
		{
			return prptProteinSeletedData;
		}
	}

	return NULL;

}

ProRataPeptideTuple * ProRataXmlProcessors::getPeptide( const QString & qsLocus,
		const QString & qsChromatogramId )
{
	setProtein( qsLocus );

	setPeptide( qsLocus, qsChromatogramId );

	if ( ! qdeCurrentPeptideElement.isNull() )
	{
		if ( processPeptide( qdeCurrentPeptideElement ) )
		{
			return prdtPeptideSeletedData;
		}
	}
	return NULL;

}

QString ProRataXmlProcessors::getText( QDomElement &qdeTag, const QString & qsTagName )
{
	return qdeTag.firstChildElement( qsTagName ).text();
}

ProRataProteinTuple * ProRataXmlProcessors::getNextProtein()
{
	if ( processProtein( qdeCurrentProteinElement ) )
	{
		qdeCurrentProteinElement = qdeCurrentProteinElement.nextSiblingElement();
		return prptProteinSeletedData;
	}
	else
	{
		return NULL;
	}

}

ProRataPeptideTuple * ProRataXmlProcessors::getNextPeptide()
{

	if ( processPeptide( qdeCurrentPeptideElement ) )
	{
		qdeCurrentPeptideElement = qdeCurrentPeptideElement.nextSiblingElement();
		return prdtPeptideSeletedData;
	}
	else
	{
		return NULL;
	}

}

void ProRataXmlProcessors::setProtein( const QString &qsLocus )
{
	if ( getText( qdeCurrentProteinElement, QString( "LOCUS" ) ).trimmed() != qsLocus.trimmed() )
	{
		qdeCurrentProteinElement = qdeFirstProteinElement;
		goToProteinNode( qsLocus );
	}

	qdeFirstPeptideElement = qdeCurrentProteinElement.firstChildElement( "CHROMATOGRAM" );
	qdeCurrentPeptideElement = qdeFirstPeptideElement;

}

void ProRataXmlProcessors::setPeptide( const QString &qsLocus, const QString &qsChromatogramId )
{
	setProtein( qsLocus );

	if ( getText( qdeCurrentPeptideElement, 
			QString( "CHROMATOGRAM_IDENTIFIER" ) ).trimmed() != 
			qsChromatogramId.trimmed() )
	{

		qdeCurrentPeptideElement = qdeCurrentProteinElement.firstChildElement( "CHROMATOGRAM" );
		goToPeptideNode( qsChromatogramId );
	}
}

void ProRataXmlProcessors::resetPointers()
{
	qdeFirstProteinElement = qdeRoot.firstChildElement( "PROTEIN" );
	qdeFirstPeptideElement = qdeFirstProteinElement.firstChildElement( "CHROMATOGRAM" );

	qdeCurrentProteinElement = qdeFirstProteinElement;
	qdeCurrentPeptideElement = qdeFirstPeptideElement;
}

bool ProRataXmlProcessors::processFile( const QString & qsFilename )
{

	QString qsErrorStr;
	int iErrorLine;
	int iErrorColumn;


	QFile qfFile( qsFilename );

	if ( !qfFile.open( QFile::ReadOnly | QFile::Text ) ) 
	{
		/*
		QMessageBox::warning( NULL, tr( "ProRata" ),
				tr( "Cannot read results file %1:\n%2." )
				.arg( fileName )
				.arg( file.errorString() ) );
		*/

		cout << "Problem with results file" << endl;

		return false; 
	}


	if ( !qddResults.setContent( &qfFile, true, &qsErrorStr, &iErrorLine,
				&iErrorColumn ) ) 
	{

		/*

		QMessageBox::information( NULL, tr("ProRata"),
				tr( "Parse error at line %1, column %2:\n%3" )
				.arg( iErrorLine )
				.arg( iErrorColumn )
				.arg( qsErrorStr ) );
		*/

		cout << "Problem parsing the results file" << endl;

		return false;
	}

	qdeRoot = qddResults.documentElement();


	if ( qdeRoot.tagName() != "QUANTITATIVE_PROTEOMICS" ||
			qdeRoot.hasAttribute( "ProRata" ) ) 
	{

		/*
		QMessageBox::information( this, tr( "ProRata" ),
				tr( "Xml file does not have any proteomics data." ) );
		*/

		cout << "Xml fiel does not have any proteomics data." << endl;

		return false;
	} 

	qdeFirstProteinElement = qdeRoot.firstChildElement( "PROTEIN" );
	qdeFirstPeptideElement = qdeFirstProteinElement.firstChildElement( "CHROMATOGRAM" );

	return true;

}

void ProRataXmlProcessors::goToProteinNode( const QString &qsLocus )
{
	while( ! qdeCurrentProteinElement.isNull() )
	{
		/*

		cout << "\"" << getText( qdeCurrentProteinElement, QString( "LOCUS" ) ).trimmed().toAscii().data() << "\" and \""
			<< qsLocus.trimmed().toAscii().data() << "\"" << endl;
		*/
		if ( getText( qdeCurrentProteinElement, QString( "LOCUS" ) ).trimmed() == qsLocus.trimmed() )
		{
			return;
		}
		qdeCurrentProteinElement = qdeCurrentProteinElement.nextSiblingElement();
	}

}

void ProRataXmlProcessors::goToPeptideNode( const QString & qsChromatogramId )
{
	while( ! qdeCurrentPeptideElement.isNull() )
	{
		if ( getText( qdeCurrentPeptideElement, QString( "CHROMATOGRAM_IDENTIFIER" ) ).trimmed() == qsChromatogramId.trimmed() )
		{
			return;
		}
		qdeCurrentPeptideElement = qdeCurrentPeptideElement.nextSiblingElement();
	}
}

