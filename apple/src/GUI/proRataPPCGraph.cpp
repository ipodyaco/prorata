
#include "proRataPPCGraph.h"

ProRataPPCGraph::ProRataPPCGraph( QWidget * qwPane ) : ProRataGraph( qwPane )
{
	// Set the title of the Graph.
	setTitle( "Paired Parallel Covariance Chromatogram" );
	
	// Set X axis label.
	setAxisTitle(QwtPlot::xBottom, "Retention time (min)");
	setAxisLabelFormat( QwtPlot::yLeft, 'e', 1, 2 );
	setAxisTitle(QwtPlot::yLeft, "Covariance");
}

ProRataPPCGraph::~ProRataPPCGraph()
{
}

void ProRataPPCGraph::cleanGraph()
{
	emit updated( QString( "PPC Chromatogram" ), 0 );
	ProRataGraph::cleanGraph();
}

void ProRataPPCGraph::setStartBoundaryMarker( double dValue )
{

	// Create a marker to show start boundary
	qwtpmStartBoundary = new QwtPlotMarker();
	qwtpmStartBoundary->setLinePen( QPen( QColor( Qt::darkGreen ), 0, Qt::DashDotLine) );
	qwtpmStartBoundary->setLabelAlignment(Qt::AlignLeft|Qt::AlignTop);
	qwtpmStartBoundary->setLineStyle(QwtPlotMarker::VLine);
	//qwtpmStartBoundary->setVisible( false );
	qwtpmStartBoundary->attach( this );
	//qwtpmStartBoundary->setLabel( ( QString( "REF @ " ) + QString::number( dValue ) ) );
	qwtpmStartBoundary->setLabel( QString::number( dValue ) );
	qwtpmStartBoundary->setLabelPen( QPen( QColor( Qt::darkGreen ) ) );
	qwtpmStartBoundary->setXValue( dValue );
	//qwtpmStartBoundary->setVisible( true );

}

void ProRataPPCGraph::setEndBoundaryMarker( double dValue )
{
	// Create a marker to show end boundary
	qwtpmEndBoundary = new QwtPlotMarker();
	qwtpmEndBoundary->setLinePen( QPen( QColor( Qt::darkGreen ) , 0, Qt::DashDotLine) );
	qwtpmEndBoundary->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmEndBoundary->setLineStyle(QwtPlotMarker::VLine);
	//qwtpmEndBoundary->setVisible( false );
	qwtpmEndBoundary->attach( this );

	//qwtpmEndBoundary->setLabel( ( QString( "TRMT @ " ) + QString::number( dValue ) ) );
	qwtpmEndBoundary->setLabel(  QString::number( dValue ) );
	qwtpmEndBoundary->setLabelPen( QPen( QColor( Qt::darkGreen ) ) );
	qwtpmEndBoundary->setXValue( dValue );
	//qwtpmEndBoundary->setVisible( true );
}

void ProRataPPCGraph::setMSMSScanMarker( double dValue )
{
	QwtPlotMarker *qwtpmPPCChro = new QwtPlotMarker;
	// QColor( 192, 131, 0 ) is brown
	qwtpmPPCChro->setLinePen( QPen( QColor( 192, 131, 0 ) , 0, Qt::DashDotLine) );
	qwtpmPPCChro->setLabelAlignment(Qt::AlignRight|Qt::AlignBottom);
	qwtpmPPCChro->setLineStyle(QwtPlotMarker::VLine);
	qwtpmPPCChro->setLabel( QString::number( dValue ) );
	qwtpmPPCChro->setLabelPen( QPen( QColor( 192, 131, 0 ) ) );
	qwtpmPPCChro->setXValue( dValue );
	qwtpmPPCChro->attach( this );
	replot();
}

void ProRataPPCGraph::setPPCData( const QwtArray<double> & dXValues, 
		const QwtArray<double> & dYValues )
{

	// Create the data interface
	qwtdPPCData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
	qwtpcPPCCurve = new QwtPlotCurve( "PPC" );
	qwtpcPPCCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcPPCCurve->setPen(QPen(Qt::black));

	// Set the data to the curve.
	qwtpcPPCCurve->setData( (*qwtdPPCData ) );

	// Add the curve to the graph
	qwtpcPPCCurve->attach( this );

	// Update the graph.
	replot();
}

void ProRataPPCGraph::peptideUpdated( PeptideRatio *peppRatio )
{

	if ( !( peppRatio ) )
	{
		QwtPlot::clear();
		replot();
		return;
	}

	QwtPlot::clear();

	vector<double> vdTime;
	vector<float> vfTime = peppRatio->getRetentionTime();

	for( unsigned int i = 0; i < vfTime.size(); i++ )
	{
		vdTime.push_back( (double)(vfTime.at(i)) );
	}

	vector<double> vdCovariannce = peppRatio->getCovarianceChro();
	double dMaxCovariance = *max_element( vdCovariannce.begin(), vdCovariannce.end() );
	double dMinCovariance = 0.0;
	setAxisScale( QwtPlot::yLeft, dMinCovariance, dMaxCovariance * 1.1 );

	setPPCData( QVector<double>::fromStdVector( vdTime ), QVector<double>::fromStdVector( vdCovariannce ) );
	
	if( fabs( peppRatio->getRightValleyTime() - peppRatio->getLeftValleyTime() ) > 0.00005 )
	{
		setStartBoundaryMarker( (double)peppRatio->getLeftValleyTime() );
		setEndBoundaryMarker( (double)peppRatio->getRightValleyTime() );
	}

	vector<float>::const_iterator vfcItr;

	for ( vfcItr = peppRatio->getMS2Time().begin();
			vfcItr != peppRatio->getMS2Time().end();
			vfcItr++ )
	{
		setMSMSScanMarker( (double)(*vfcItr) );
	}

	emit updated( QString( "PPC Chromatogram" ), 1 );
	replot();
}
