
#include "proRataLikelihoodGraph.h"
#include <QtGlobal>
#include <iostream>

using namespace std;

ProRataLikelihoodGraph::ProRataLikelihoodGraph( QWidget * qwPane ) : ProRataGraph( qwPane )
{

	// Set the title of the Graph.
	setTitle( "Profile Likelihood Curve" );

	// Disable all axes.
	enableYRightAxis( true );
	enableYLeftAxis( true );

	setAxisAutoScale(QwtPlot::yLeft);
	setAxisAutoScale( QwtPlot::yRight );
	
//	qwtppPicker->setSelectionFlags(QwtPicker::RectSelection | QwtPicker::DragSelection);	
//	qwtppPicker->setRubberBand(QwtPicker::RectRubberBand);
	
	qwtppPicker->setSelectionFlags(QwtPicker::PointSelection | QwtPicker::ClickSelection);	
	
//	connect( qwtppPicker, SIGNAL( selected( const QwtDoubleRect & ) ),
//			this, SLOT( relayRectSelected( const QwtDoubleRect & ) ) );
	connect( qwtppPicker, SIGNAL( selected ( const QwtDoublePoint & ) ),
			this, SLOT( relayPointSelected( const QwtDoublePoint & ) ) );	// Set X axis label.
	setAxisTitle(QwtPlot::xBottom, "Log2 (Ratio)");

	bPreviousSelected = false;

	bPopulated = false;

	qwtpcSelectedPeptide = NULL;


}

ProRataLikelihoodGraph::~ProRataLikelihoodGraph()
{
}

void ProRataLikelihoodGraph::cleanGraph()
{
	emit updated( QString( "Likelihood Plot" ), 0 );
	ProRataGraph::cleanGraph();
}

void ProRataLikelihoodGraph::setLowerCI( double dValue )
{
	// Create a marker to show lower confidence interval.
	qwtpmLowerCIMarker = new QwtPlotMarker();
	qwtpmLowerCIMarker->setLinePen( QPen( Qt::black , 0, Qt::DashDotLine) );
	qwtpmLowerCIMarker->setLabelAlignment(Qt::AlignLeft|Qt::AlignTop);
	qwtpmLowerCIMarker->setLineStyle(QwtPlotMarker::VLine);

	qwtpmLowerCIMarker->setLabel(  QString::number( dValue )  );
	qwtpmLowerCIMarker->setXValue( dValue );
	qwtpmLowerCIMarker->attach( this );

}

void ProRataLikelihoodGraph::setUpperCI( double dValue )
{
	// Create a marker to show upper confidence interval.
	qwtpmUpperCIMarker = new QwtPlotMarker();
	qwtpmUpperCIMarker->setLinePen( QPen( Qt::black , 0, Qt::DashDotLine) );
	qwtpmUpperCIMarker->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmUpperCIMarker->setLineStyle(QwtPlotMarker::VLine);
	//qwtpmUpperCIMarker->setVisible( false );

	qwtpmUpperCIMarker->setLabel(  QString::number( dValue )  );
	qwtpmUpperCIMarker->setXValue( dValue );
	qwtpmUpperCIMarker->attach( this );
	//qwtpmUpperCIMarker->setVisible( true );
}

void ProRataLikelihoodGraph::setCutoff( double dCutoff )
{

	// Create a horizontal line to show the cutoff value.
	qwtpmLikelihoodMarker = new QwtPlotMarker();
	qwtpmLikelihoodMarker->setLinePen( QPen( Qt::black , 0, Qt::DashDotLine) );
	qwtpmLikelihoodMarker->setLabelAlignment(Qt::AlignRight|Qt::AlignBottom);
	qwtpmLikelihoodMarker->setLineStyle(QwtPlotMarker::HLine);
	qwtpmLikelihoodMarker->setYAxis(QwtPlot::yRight);

	qwtpmLikelihoodMarker->setLabel(  QString::number( dCutoff )  );
	qwtpmLikelihoodMarker->setYValue( dCutoff );
	qwtpmLikelihoodMarker->attach( this );

}

void ProRataLikelihoodGraph::setLikelihoodData( const QwtArray<double> & dXValues, 
		const QwtArray<double> &dYValues )
{
	// Update Y axis label
	setAxisTitle(QwtPlot::yRight, "Ln (Likelihood)");
	setAxisLabelFormat( QwtPlot::yRight, 'f', 1, 3 );
	// Create the data interface
	//qwtdLikelihoodData = new QwtArrayData( dXValues, dYValues, iSize );
	qwtdLikelihoodData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
	qwtpcLikelihoodCurve = new QwtPlotCurve( "Likelihood Curve" );
	qwtpcLikelihoodCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcLikelihoodCurve->setStyle( QwtPlotCurve::Lines );
	qwtpcLikelihoodCurve->setYAxis(QwtPlot::yRight);
	QPen qpPen(Qt::blue);
	qpPen.setWidth( 1 );
	qpPen.setJoinStyle( Qt::RoundJoin );
	qwtpcLikelihoodCurve->setPen( qpPen );
//	qwtpcLikelihoodCurve->setPen(QPen(Qt::blue));

	// Set the data to the curve.
	qwtpcLikelihoodCurve->setData( (*qwtdLikelihoodData ) );

	// Add the curve to the graph
	qwtpcLikelihoodCurve->attach( this );

	
	// Update the graph.
	replot();

}

void ProRataLikelihoodGraph::setPeptidePoints( const QwtArray<double> & dXValues, 
		const QwtArray<double> &dYValues )
{

	/*
	// Enable Axis
	enableYLeftAxis( true );
	 */

	// Update Y axis label
	setAxisTitle(QwtPlot::yLeft, "Log2 (SNR)" );
	setAxisLabelFormat( QwtPlot::yLeft, 'f', 1, 2 );

	// Create the data inteface.
	//qwtdPeptidePointsData = new QwtArrayData( dXValues, dYValues, iSize );
	qwtdPeptidePointsData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the data.
	qwtpcPeptidePointsCurve = new QwtPlotCurve( "Peptide data points" );

	// Create a symbol to represent data  points.
	qwtsPeptidePointsSym.setStyle(QwtSymbol::Ellipse);
	qwtsPeptidePointsSym.setPen(QColor(Qt::red));
	qwtsPeptidePointsSym.setBrush(QBrush(Qt::red));
	qwtsPeptidePointsSym.setSize(4);

	qwtpcPeptidePointsCurve->setSymbol( qwtsPeptidePointsSym );
	qwtpcPeptidePointsCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcPeptidePointsCurve->setYAxis(QwtPlot::yLeft);
	// This curve will not have any connection between points.
	qwtpcPeptidePointsCurve->setStyle( QwtPlotCurve::NoCurve );

	// Add the data to the curve.
	qwtpcPeptidePointsCurve->setData( (*qwtdPeptidePointsData ) );

	// Add the curve to the graph
	qwtpcPeptidePointsCurve->attach( this );

	
	// Update the graph.
	replot();

}

void ProRataLikelihoodGraph::proteinUpdated( ProteinRatio *propRatio )
{
	if (!(propRatio) )
	{
		QwtPlot::clear();
		replot();
		return;
	}

	propRatioCurrent = propRatio;

	QwtPlot::clear();
	bPreviousSelected = false;
	
	setLikelihoodData( 
		QVector<double>::fromStdVector( propRatio->getLog2RatioBin() ),
		QVector<double>::fromStdVector( propRatio->getLnLikelihood() ) );

	setPeptidePoints( 
		QVector<double>::fromStdVector( propRatio->getPeptideLog2Ratio() ),
		QVector<double>::fromStdVector( propRatio->getPeptideLog2SN() ) );

	setLowerCI( propRatio->getLowerLimitCI() );
	setUpperCI( propRatio->getUpperLimitCI() );
	setCutoff( propRatio->getLnLikelihoodCutoff() );

	/*emit updated( LIKELIHOOD );*/
	bPopulated = true;
	emit updated( QString( "Likelihood Plot" ), 1 );
	replot();
}
/*
void ProRataLikelihoodGraph::relayRectSelected( const QwtDoubleRect &rect )
{
	emit rectMLESelected( rect );
}
*/

void ProRataLikelihoodGraph::relayPointSelected( const QwtDoublePoint & point )
{
	if (!bPopulated)
	{
		return;
	}
	emit pointMLESelected( point );
}

void ProRataLikelihoodGraph::peptideUpdated( PeptideInfo * pPepInfo )
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


	double dPeptideLogSN = pPepInfo->getPCALog2SNR();
	double dPeptideLogRatio = pPepInfo->getPCALog2Ratio();
	
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
	
	qwtpcSelectedPeptide = new QwtPlotCurve( "Selected Peptide" );
	bPreviousSelected = true;
	
//	qDebug( "inside slot %f", dPeptideLogSN );

	// Create a symbol to represent data  points.
	qwtsPeptideSelectedSym.setStyle(QwtSymbol::Ellipse);
	qwtsPeptideSelectedSym.setPen(QColor(Qt::green));
	qwtsPeptideSelectedSym.setBrush(QBrush(Qt::green));
	qwtsPeptideSelectedSym.setSize(6);

	qwtpcSelectedPeptide->setSymbol( qwtsPeptideSelectedSym );
	qwtpcSelectedPeptide->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcSelectedPeptide->setYAxis(QwtPlot::yLeft);
	// This curve will not have any connection between points.
	qwtpcSelectedPeptide->setStyle( QwtPlotCurve::NoCurve );

	// Add the data to the curve.
	qwtpcSelectedPeptide->setData( &dPeptideLogRatio, &dPeptideLogSN, 1 );

	// Add the curve to the graph
	qwtpcSelectedPeptide->attach( this );

	// Update the graph.
	replot();

}
