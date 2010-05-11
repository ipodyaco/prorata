
#include "proRataTextArea.h"

#include <iostream>
#include <QMessageBox>
#include <math.h>
#include <QtGlobal>

using namespace std;

ProRataTextArea::ProRataTextArea( QWidget * qwParent ) :
	QTextBrowser( qwParent )
{
	QString qsHTMLpepCurrent = "";
	QString qsHTMLproCurrent = "";

	QString qsHTMLpepOrigninal = "";
	QString qsHTMLproOrigninal = "";

	QString qsHTMLContentCurrent = "";
	QString qsHTMLContentOriginal = "";	
	loadContent();

}

ProRataTextArea::~ProRataTextArea()
{
}

void ProRataTextArea::loadContent()
{
	qsHTMLproOrigninal = QString( "<html><body><h4>Protein Information </h4><table><tr><td bgcolor=#dddbdb>Protein description</td><td><!--Prop desc-->--Desciption--<!--Prop desc--></td></tr><tr><td bgcolor=#dddbdb>Ratio of samples</td><td><!--tre to ref-->--Numerator : Denominator--<!--tre to ref--></td></tr></table><table><tr><td bgcolor=#dddbdb>Abundance ratio</td><td><!--Pro abundance ratio-->--Ratio--<!--Pro abundance ratio--></td><td bgcolor=#dddbdb>Log2(abundance ratio)</td><td><!--Pro log2 ratio-->--Log2 Ratio--<!--Pro log2 ratio--></td><td bgcolor=#dddbdb>Confidence interval</td><td><!--Pro CI-->--[Minimum, Maximum]--<!--Pro CI--></td></tr><tr><td bgcolor=#dddbdb>Quantification hits</td><td><!--Quant hit count-->--Count--<!--Quant hit count--></td><td bgcolor=#dddbdb>Identification hits</td><td><!--ID hit count-->--Count--<!--ID hit count--></td><td bgcolor=#dddbdb>Locus</td><td><!--prot locus-->--Protein--<!--prot locus--></td></tr></table>" );
	qsHTMLpepOrigninal = QString( "<h4>Peptide Information</h4><table><tr><td bgcolor=#dddbdb>Sequence</td><td><!--Sequence-->--Sequence--<!--Sequence--></td></tr><tr><td bgcolor=#dddbdb>Identification filename</td><td><!--ID fn-->--File Name--<!--ID fn--></td></tr><tr><td bgcolor=#dddbdb>Chromatogram filename</td><td><!--Chro fn-->--File Name--<!--Chro fn--></td></tr></table><table><tr><td bgcolor=#dddbdb>XIC identifier</td><td><!--Chro number-->--Number--<!--Chro number--></td><td bgcolor=#dddbdb>Peak shift</td><td><!--MS1 shifted-->--Shifted Scans--<!--MS1 shifted--></td><td bgcolor=#dddbdb>ID score</td><td><!--max id score-->--Maximum Score--<!--max id score--></td></tr><tr><td bgcolor=#dddbdb>PCA abundance ratio</td><td><!--PCA ratio-->--Abundance Ratio--<!--PCA ratio--></td><td bgcolor=#dddbdb>Peak area ratio</td><td><!--area ratio-->--Abundance Ratio--<!--area ratio--></td><td bgcolor=#dddbdb>Peak height ratio</td><td><!--height ratio-->--Abundance Ratio--<!--height ratio--></td></tr><tr><td bgcolor=#dddbdb>Profile SNR</td><td><!--EV SNR-->--S/N--<!--EV SNR--></td><td bgcolor=#dddbdb>Reference chro SNR</td><td><!--ref SNR-->--S/N--<!--ref SNR--></td><td bgcolor=#dddbdb>Treatment chro SNR</td><td><!--tre SNR-->--S/N--<!--tre SNR--></td></tr><tr><td bgcolor=#dddbdb>Full Scan Range</td><td><!--MS1 range-->--scan range--<!--MS1 range--></td><td bgcolor=#dddbdb>Treatment m/z range</td><td><!--tre mz-->--m/z range--<!--tre mz--></td><td bgcolor=#dddbdb>Reference m/z range</td><td><!--ref mz-->--m/z range--<!--ref mz--></td></tr></table></body></html>" );
	qsHTMLproOrigninal.replace(QString("Treatment"), QString( ProRataConfig::getNumeratorIsotopologue().c_str() )  );
	qsHTMLproOrigninal.replace(QString("Reference"), QString( ProRataConfig::getDenominatorIsotopologue().c_str() ) );
	qsHTMLpepOrigninal.replace(QString("Treatment"), QString( ProRataConfig::getNumeratorIsotopologue().c_str() )  );
	qsHTMLpepOrigninal.replace(QString("Reference"), QString( ProRataConfig::getDenominatorIsotopologue().c_str() ) );
	QString qsRatioName =  QString( ProRataConfig::getNumeratorIsotopologue().c_str() )
	       	+ QString( " : " ) + QString( ProRataConfig::getDenominatorIsotopologue().c_str() );
	substituteContent(qsHTMLproOrigninal, QString( "<!--tre to ref-->" ), qsRatioName );
	qsHTMLpepCurrent = qsHTMLpepOrigninal;
	qsHTMLproCurrent = qsHTMLproOrigninal;
	qsHTMLContentCurrent = qsHTMLproOrigninal + qsHTMLpepOrigninal;
	insertHtml( qsHTMLContentCurrent );
}

void ProRataTextArea::substituteContent( QString & qsHTMLContent, const QString & qsTag, const QString & qsWithValue )
{
	QString qsTemp = qsHTMLContent;

	QString qsReplacementString = qsTemp.section( qsTag, 1, 1 );
	qsReplacementString = qsTag + qsReplacementString + qsTag;

	QString qsNewString = qsTag + qsWithValue + qsTag;

	qsHTMLContent.replace( qsReplacementString, qsNewString );
}

void ProRataTextArea::updateProteinDesc( ProteinInfo * proInfo )
{
	if (!(proInfo) )
	{
		qsHTMLContentCurrent = qsHTMLproOrigninal + qsHTMLpepOrigninal;
		clear();
		insertHtml( qsHTMLContentCurrent );
		return;
	}

	// update protein description
	substituteContent( qsHTMLproCurrent, QString( "<!--Prop desc-->" ), QString( proInfo->getDescription().c_str() ) );

	// update protein locus
	substituteContent( qsHTMLproCurrent, QString( "<!--prot locus-->" ), QString( proInfo->getLocus().c_str() ) );
	
	// update confidence interval
	QString qsCI = "[";
	qsCI.append( QString::number(proInfo->getLowerLimitCI() ) );
	qsCI.append(  ", " );
	qsCI.append( QString::number(proInfo->getUpperLimitCI() ) );
	qsCI.append( "]" );
	substituteContent( qsHTMLproCurrent, QString( "<!--Pro CI-->" ), qsCI );

	// update abundance ratio
	substituteContent( qsHTMLproCurrent, QString( "<!--Pro log2 ratio-->" ), QString::number( proInfo->getLog2Ratio() ) );
	substituteContent( qsHTMLproCurrent, QString( "<!--Pro abundance ratio-->" ), QString::number( pow( 2.0, proInfo->getLog2Ratio() ), 'g', 2 ) );

	// update the hit count
	substituteContent( qsHTMLproCurrent, QString( "<!--Quant hit count-->" ), QString::number( proInfo->getQuantifiedPeptides() ) );
	substituteContent( qsHTMLproCurrent, QString( "<!--ID hit count-->" ), QString::number( proInfo->getIdentifiedPeptides() ) );
	
	// remove the updates in the peptide content
	qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepOrigninal;
	
	clear();
	insertHtml( qsHTMLContentCurrent );
}


void ProRataTextArea::updateProteinDesc( ProteinRatio * proRatio )
{
//	qDebug( "Inside updateProteinDesc( ProteinRatio * proRatio )" );
}

void ProRataTextArea::updatePeptideDesc( PeptideInfo * pepInfo )
{
//	qDebug( "Inside updatePeptideDesc( PeptideInfo * pepInfo )" );

	if ( ! ( pepInfo ) )
	{

		qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepOrigninal;
		clear();
		insertHtml( qsHTMLContentCurrent );
		return;
	}

	substituteContent( qsHTMLpepCurrent, QString( "<!--Chro fn-->" ), QString( pepInfo->getFilename().c_str() ) );

	qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepCurrent;
	clear();
	insertHtml( qsHTMLContentCurrent );

}

void ProRataTextArea::updatePeptideDesc( PeptideRatio * pepRatio )
{

	if ( ! ( pepRatio ) )
	{

		qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepOrigninal;
		clear();
		insertHtml( qsHTMLContentCurrent );
		return;
	}

	unsigned int i;
	// update the sequence
	substituteContent( qsHTMLpepCurrent, QString( "<!--Sequence-->" ),  QString::fromStdString( pepRatio->getSequence() ) );

	// update the chromatogram identifier	
	substituteContent( qsHTMLpepCurrent, QString( "<!--Chro number-->" ),  QString::number( pepRatio->getIdentifier() ) );

	// update the peak shift	
	substituteContent( qsHTMLpepCurrent, QString( "<!--MS1 shifted-->" ),  QString::number( pepRatio->getPeakShift() ) );

	// update the ID score	
	substituteContent( qsHTMLpepCurrent, QString( "<!--max id score-->" ),  QString::number( pepRatio->getMaximumScore() ) );
	
	// update the quantification results
	substituteContent( qsHTMLpepCurrent, QString( "<!--PCA ratio-->" ),  QString::number( pepRatio->getPCARatio() ) );
	substituteContent( qsHTMLpepCurrent, QString( "<!--area ratio-->" ),  QString::number( pepRatio->getPeakAreaRatio() ) );
	substituteContent( qsHTMLpepCurrent, QString( "<!--height ratio-->" ),  QString::number( pepRatio->getPeakHeightRatio() ) );
	substituteContent( qsHTMLpepCurrent, QString( "<!--EV SNR-->" ),  QString::number( pepRatio->getPCASN() ) );
	substituteContent( qsHTMLpepCurrent, QString( "<!--ref SNR-->" ),  QString::number( pepRatio->getReferencePeakSNR() ) );
	substituteContent( qsHTMLpepCurrent, QString( "<!--tre SNR-->" ),  QString::number( pepRatio->getTreatmentPeakSNR() ) );
	
	// update the scan range
	QString qsScanRange = "[";
	qsScanRange.append( QString::number( pepRatio->getFirstScanNumber() ) );
	qsScanRange.append(  ", " );
	qsScanRange.append( QString::number( pepRatio->getLastScanNumber() ) );
	qsScanRange.append( "]" );
	substituteContent( qsHTMLpepCurrent, QString( "<!--MS1 range-->" ),  qsScanRange );
	
	// update the mz range
	QString qsReferenceMZ = "";
	QString qsTreatmentMZ = "";
	
	vector< float > vfLowerMZ;
	vector< float > vfUpperMZ;
	
	if( pepRatio->getReferenceMZrange( vfLowerMZ, vfUpperMZ ) )
	{
		for( i = 0; i < vfLowerMZ.size(); ++i )
		{
			qsReferenceMZ.append( "[" );
			qsReferenceMZ.append(  QString::number( vfLowerMZ[i] ) );
			qsReferenceMZ.append(  ", " );
			qsReferenceMZ.append(  QString::number( vfUpperMZ[i] ) );
			qsReferenceMZ.append( "]" );
			qsReferenceMZ.append( "\n" );
		}
	}
	substituteContent( qsHTMLpepCurrent, QString( "<!--ref mz-->" ), qsReferenceMZ );
	
	if( pepRatio->getTreatmentMZrange( vfLowerMZ, vfUpperMZ ) )
	{
		for( i = 0; i < vfLowerMZ.size(); ++i )
		{
			qsTreatmentMZ.append( "[" );
			qsTreatmentMZ.append(  QString::number( vfLowerMZ[i] ) );
			qsTreatmentMZ.append(  ", " );
			qsTreatmentMZ.append(  QString::number( vfUpperMZ[i] ) );
			qsTreatmentMZ.append( "]" );
			qsTreatmentMZ.append( "\n" );
		}
	}
	substituteContent( qsHTMLpepCurrent, QString( "<!--tre mz-->" ), qsTreatmentMZ );

	// update the ID filenames
	string sFilenameList = "";
	vector< string > vsFilename = pepRatio->getAllIDfilename();

	qDebug(  "size %i", vsFilename.size()  );
	for( i = 0; i < vsFilename.size() ; ++i )
	{
		qDebug(  "filename"  );
		qDebug(  vsFilename[i].c_str()  );
		sFilenameList.append( vsFilename[i] );
		sFilenameList.append( "<br>" );
	}
	substituteContent( qsHTMLpepCurrent, QString( "<!--ID fn-->" ),  QString::fromStdString(sFilenameList) );

	// update HTML content
	qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepCurrent;
	clear();
	insertHtml( qsHTMLContentCurrent );
}

void ProRataTextArea::cleanUp()
{
	qsHTMLContentCurrent = qsHTMLproOrigninal + qsHTMLpepOrigninal;
	clear();
	insertHtml( qsHTMLContentCurrent );
	return;
}

void ProRataTextArea::cleanUpPeptide()
{
	qsHTMLContentCurrent = qsHTMLproCurrent + qsHTMLpepOrigninal;
	clear();
	insertHtml( qsHTMLContentCurrent );
	return;
}


