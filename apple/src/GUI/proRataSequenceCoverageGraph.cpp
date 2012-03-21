
#include "proRataSequenceCoverageGraph.h"

#include <QtGlobal>
#include <iostream>
#include <QVector>
using namespace std;

ProRataSequenceCoverageGraph::ProRataSequenceCoverageGraph( QWidget * qwPane ) : ProRataGraph( qwPane )
{
	// Set the title of the Graph.
	setTitle( "Quantitative Sequence Coverage" );

	setAxisTitle(QwtPlot::xBottom, "Log2 (Ratio)" );
	setAxisTitle(QwtPlot::yLeft, "Residue Number");

//	legend()->setDisplayPolicy( QwtLegend::None );

	pFASTAdata = new FASTAdata;
	bool bReadFASTA = pFASTAdata->readFASTA( ProRataConfig::getFASTAfilename() );
	if( !bReadFASTA )
		qDebug( "cannot read fASTA ");

	setXAxisRange( ProRataConfig::getMLEMinLog2Ratio(), ProRataConfig::getMLEMaxLog2Ratio() );
	bPreviousSelected = false;
	qwtpcSelectedPeptide = NULL;
	sProteinSequence = "";
	
}

ProRataSequenceCoverageGraph::~ProRataSequenceCoverageGraph()
{
	delete pFASTAdata;
}

void ProRataSequenceCoverageGraph::cleanGraph()
{
	emit updated( QString( "Seq Cov" ), 0 );
	ProRataGraph::cleanGraph();
}

void ProRataSequenceCoverageGraph::proteinUpdated( ProteinInfo * pProteinInfo )
{

	if (!(pProteinInfo) )
	{
		QwtPlot::clear();
		replot();
		return;
	}

	QwtPlot::clear();
	bPreviousSelected = false;
	
	sProteinSequence = pFASTAdata->getProteinSequence( ( pProteinInfo->getLocus() ) );
//	qDebug( "after Protein sequence %d", sProteinSequence.length() );
	setYAxisRange( 1, sProteinSequence.size() );
	vector< PeptideInfo * > vPeptideInfo = pProteinInfo->getPeptideInfo();
//	qDebug( "after peptide info " );

	for( unsigned int i = 0; i < vPeptideInfo.size(); ++i )
	{
		if( vPeptideInfo[i]->getValidity() )
		{
			double dPeptidelog2R = vPeptideInfo[i]->getPCALog2Ratio();
			int iNcord = 0;
			int iCcord = 0;
			pFASTAdata->computeCoveragePlot( sProteinSequence, vPeptideInfo[i]->getSequence(),
					iNcord, iCcord );
			addPeptide( dPeptidelog2R, iNcord, iCcord );
		}
	}
	
	// Update the graph.
	emit updated( QString( "Seq Cov" ), 1 );
	replot();

}

void ProRataSequenceCoverageGraph::addPeptide( double dX, int iY1, int iY2 )
{

	QwtArray<double> qwtdXData;
	QwtArray<double> qwtdYData;

	qwtdXData.push_back( dX );
	qwtdYData.push_back( (double)iY1 );
	qwtdXData.push_back( dX );
	qwtdYData.push_back( (double)iY2 );

	// Create the data interface
	QwtArrayData qwtdNewData( qwtdXData, qwtdYData );

	// Create a curve to represent the given data.
	

	QwtPlotCurve *qwtpcPeptide = new QwtPlotCurve;
	qwtpcPeptide->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcPeptide->setStyle( QwtPlotCurve::Lines );
	qwtpcPeptide->setPen(QPen(Qt::red));


	// Set the data to the curve.
	qwtpcPeptide->setData( qwtdNewData );

	// Add the curve to the graph
	qwtpcPeptide->attach( this );

//	legend()->hide();

//	legendStatus( false );



	//legendStatus( false );
}

void ProRataSequenceCoverageGraph::setXAxisRange( double dLow, double dHigh )
{
	setAxisScale( xBottom, dLow, dHigh );
}

void ProRataSequenceCoverageGraph::setYAxisRange( int iLow, int iHigh )
{
	setAxisScale( yLeft, (double)iLow, (double)iHigh );
}


void ProRataSequenceCoverageGraph::peptideUpdated( PeptideInfo * pPepInfo )
{
	if ( !( pPepInfo ) )
	{
		if ( qwtpcSelectedPeptide )
		{
			qwtpcSelectedPeptide->hide();
		}
		replot();
		bPreviousSelected = false;
		return;
	}
	
	qDebug( "after peptide select %d", sProteinSequence.length() );
	
	if( bPreviousSelected )
	{
		delete qwtpcSelectedPeptide;
		qwtpcSelectedPeptide = NULL;
	}

	if( !pPepInfo->getValidity() )
	{
		replot();
		bPreviousSelected = false;
		return;
	}

	qDebug( "after validity %d", sProteinSequence.length() );
	
	double dPeptidelog2R = pPepInfo->getPCALog2Ratio();
	int iNcord = 0;
	int iCcord = 0;
	pFASTAdata->computeCoveragePlot( sProteinSequence, pPepInfo->getSequence(),
			iNcord, iCcord );
	qDebug( "after FASTAy %d", iNcord );
//	
//	use the code of the addPeptide function to add the curve.
//	
	double dX = dPeptidelog2R;
	int iY1 = iNcord;
	int iY2 = iCcord;
	
	QwtArray<double> qwtdXData;
	QwtArray<double> qwtdYData;

	qwtdXData.push_back( dX );
	qwtdYData.push_back( (double)iY1 );
	qwtdXData.push_back( dX );
	qwtdYData.push_back( (double)iY2 );

	// Create the data interface
	QwtArrayData qwtdNewData( qwtdXData, qwtdYData );

	// Create a curve to represent the given data.
	
	
	qwtpcSelectedPeptide = new QwtPlotCurve( "Selected Peptide" );
	bPreviousSelected = true;
	
	qDebug( "after creat %d", iNcord );

	qwtpcSelectedPeptide->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcSelectedPeptide->setStyle( QwtPlotCurve::Lines );
	QPen qpPen(Qt::green);
	qpPen.setWidth( 2 );
	qwtpcSelectedPeptide->setPen(qpPen);


	// Set the data to the curve.
	qwtpcSelectedPeptide->setData( qwtdNewData );

	// Add the curve to the graph
	qwtpcSelectedPeptide->attach( this );

	qDebug( "after attacht %d", iNcord );
	// Update the graph.
	replot();

}

