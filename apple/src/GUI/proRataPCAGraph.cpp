
#include "proRataPCAGraph.h"
#include <QtGlobal>
#include "qwt_plot_layout.h"

ProRataPCAGraph::ProRataPCAGraph( QWidget * qwPane ) : ProRataGraph( qwPane )
{

	qwtdPeptideData = NULL;
	
	// Set the title of the Graph.
	setTitle( "Principle Component Analysis" );

	// Set axis label.
	setAxisTitle(QwtPlot::xBottom,  QString( "Intensity " ) + QString( ProRataConfig::getDenominatorIsotopologue().c_str() ) );
	setAxisTitle(QwtPlot::yLeft,  QString( "Intensity " ) + QString( ProRataConfig::getNumeratorIsotopologue().c_str() ) );

	setAxisLabelFormat( QwtPlot::xBottom, 'e', 1, 2 );
	setAxisLabelFormat( QwtPlot::yLeft, 'e', 1, 2 );

	peppRatioCurrent = NULL;
	qwtpcSelectedMS1 = NULL;

	bSecondPCA = false;
	bIsMS1Selected = false;

}

ProRataPCAGraph::~ProRataPCAGraph()
{
}

void ProRataPCAGraph::cleanGraph()
{
	emit updated( QString( "PCA Plot" ), 0 );
	ProRataGraph::cleanGraph();
}

void ProRataPCAGraph::setPeptideData( const QwtArray<double> & dXValues, 
		const QwtArray<double> & dYValues)
{
	// Create the data interface
	qwtdPeptideData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
	qwtpcPeptideCurve = new QwtPlotCurve( "MS full scans" );
	qwtpcPeptideCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcPeptideCurve->setStyle( QwtPlotCurve::NoCurve );
	qwtpcPeptideCurve->setPen(QPen(Qt::blue));

	// Create a symbol to represent data  points.
	qwtsPeptideSym.setStyle(QwtSymbol::Ellipse);
	qwtsPeptideSym.setPen(QColor(Qt::blue));
	qwtsPeptideSym.setBrush(QBrush(Qt::blue));
	qwtsPeptideSym.setSize(4);

	qwtpcPeptideCurve->setSymbol( qwtsPeptideSym );

	// Set the data to the curve.
	qwtpcPeptideCurve->setData( (*qwtdPeptideData ) );

	// Add the curve to the graph
	qwtpcPeptideCurve->attach( this );

}


void ProRataPCAGraph::setPrincipalComponent( double dX1, double dY1,
				double dX2, double dY2 )
{

	QwtArray<double> qwtdXData;
	QwtArray<double> qwtdYData;

	qwtdXData.push_back( dX1 );
	qwtdXData.push_back( dX2 );
	qwtdYData.push_back( dY1 );
	qwtdYData.push_back( dY2 );

	QwtArrayData *qwtdNewData = new QwtArrayData( qwtdXData, qwtdYData );

	QwtPlotCurve *qwtpcFirstPC;

	QPen qpPen(Qt::blue);
	qpPen.setWidth( 2 );
	if ( bSecondPCA)
	{
		qpPen.setColor(Qt::black);
		qwtpcFirstPC = new QwtPlotCurve( "PC1" );
		qwtpcFirstPC->setPen( qpPen );
	}
	else
	{
		qpPen.setColor(Qt::red);
		qwtpcFirstPC = new QwtPlotCurve( "PC2" );
		qwtpcFirstPC->setPen(qpPen);
	}
	
	qwtpcFirstPC->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcFirstPC->setStyle( QwtPlotCurve::Lines );
	

	// Set the data to the curve.
	qwtpcFirstPC->setData( (*qwtdNewData ) );

	// Add the curve to the graph
	qwtpcFirstPC->attach( this );

}

bool ProRataPCAGraph::rangeOfData( double *dYMin, double *dYMax, double *dXMin, 
		double *dXMax, double *dXAvg, double *dYAvg )
{
	if ( !qwtdPeptideData )
	{
		return false;
	}

	if ( qwtdPeptideData->size() <= 0 )
	{
		return false;
	}

	*dYMin = qwtdPeptideData->y(0);
	*dYMax = qwtdPeptideData->y(0);

	*dXMin = qwtdPeptideData->x(0);
	*dXMax = qwtdPeptideData->x(0);

	*dXAvg = qwtdPeptideData->x(0);
	*dYAvg = qwtdPeptideData->y(0);

	for ( int i = 1; i < (int)qwtdPeptideData->size(); i++ )
	{
		if ( (*dYMin) < qwtdPeptideData->y(i) )
		{	
			*dYMin = qwtdPeptideData->y(i);
		}

		if ( (*dYMax) > qwtdPeptideData->y(i) )
		{	
			*dYMax = qwtdPeptideData->y(i);
		}

		if ( (*dXMin) < qwtdPeptideData->x(i) )
		{	
			*dXMin = qwtdPeptideData->x(i);
		}

		if ( (*dXMax) > qwtdPeptideData->x(i) )
		{	
			*dXMax = qwtdPeptideData->x(i);
		}

		*dXAvg = *dXAvg + qwtdPeptideData->x(i);
		*dYAvg = *dYAvg + qwtdPeptideData->y(i);
	}

	*dXAvg = (*dXAvg) / qwtdPeptideData->size();
	*dYAvg = (*dYAvg) / qwtdPeptideData->size();

	return true;
}

void ProRataPCAGraph::peptideUpdated( PeptideRatio *peppRatio )
{

	if ( !( peppRatio ) )
	{
		QwtPlot::clear();
		bIsMS1Selected = false;
		replot();
		return;
	}

	peppRatioCurrent = peppRatio;

	bIsMS1Selected = false;
	QwtPlot::clear();
	
	vector< double > vdDataX = peppRatio->getReferencePeak();
	vector< double > vdDataY = peppRatio->getTreatmentPeak();

	setPeptideData( 
		QVector<double>::fromStdVector( vdDataX ),
		QVector<double>::fromStdVector( vdDataY ) );
	double dXAve, dYAve, dYMin, dYMax, dXMin, dXMax;
	if ( !rangeOfData( &dYMin, &dYMax, &dXMin, &dXMax, &dXAve, &dYAve ) )
		return;

	vector< double > vdAxisPoints;
	vdAxisPoints.push_back( dYMin );
	vdAxisPoints.push_back( dYMax );
	vdAxisPoints.push_back( dXMin );
	vdAxisPoints.push_back( dXMax );
	
	
	double dPC1x, dPC1y, dPC2x, dPC2y;
	peppRatio->getEngenvectors( dPC1x, dPC1y, dPC2x, dPC2y );
		
	double dPC1EV, dPC2EV;
	peppRatio->getEngenvalues( dPC1EV, dPC2EV );


	double dLength = sqrt( (dXMax - dXAve)*(dXMax - dXAve) + (dYMax - dYAve)*(dYMax - dYAve) );
	
	double dPC1Length = dLength * dPC1EV / ( dPC1EV + dPC2EV) ;
	double dPC1LengthUnit = sqrt( dPC1x*dPC1x + dPC1y*dPC1y );
	double dPC1x1 =  (dPC1x/dPC1LengthUnit) * dPC1Length + dXAve;
	double dPC1y1 =  (dPC1y/dPC1LengthUnit) * dPC1Length + dYAve;
	double dPC1x2 = -(dPC1x/dPC1LengthUnit) * dPC1Length + dXAve;
	double dPC1y2 = -(dPC1y/dPC1LengthUnit) * dPC1Length + dYAve;
	bSecondPCA = false;
	setPrincipalComponent( dPC1x1, dPC1y1, dPC1x2, dPC1y2 );
	vdAxisPoints.push_back( dPC1x1 );
	vdAxisPoints.push_back( dPC1x2 );
	vdAxisPoints.push_back( dPC1y1 );
	vdAxisPoints.push_back( dPC1y2 );

	double dPC2Length = dLength * dPC2EV / ( dPC1EV + dPC2EV);
	double dPC2LengthUnit = sqrt( dPC2x*dPC2x + dPC2y*dPC2y );
	double dPC2x1 =  (dPC2x/dPC2LengthUnit) * dPC2Length + dXAve;
	double dPC2y1 =	 (dPC2y/dPC2LengthUnit) * dPC2Length + dYAve;
	double dPC2x2 = -(dPC2x/dPC2LengthUnit) * dPC2Length + dXAve;
	double dPC2y2 = -(dPC2y/dPC2LengthUnit) * dPC2Length + dYAve;
	bSecondPCA = true;
	setPrincipalComponent( dPC2x1, dPC2y1, dPC2x2, dPC2y2 );
	vdAxisPoints.push_back( dPC2x1 );
	vdAxisPoints.push_back( dPC2x2 );
	vdAxisPoints.push_back( dPC2y1 );
	vdAxisPoints.push_back( dPC2y2 );

	double dAxisMax = *max_element( vdAxisPoints.begin(), vdAxisPoints.end() );
	double dAxisMin = *min_element( vdAxisPoints.begin(), vdAxisPoints.end() );
	setAxisScale( QwtPlot::yLeft, dAxisMin * 0.9, dAxisMax * 1.1 );
	setAxisScale( QwtPlot::xBottom, dAxisMin * 0.9, dAxisMax * 1.1 );
	
//	qDebug( "dPC1EV: %f", dPC1EV );
//	qDebug( "dPC2EV: %f", dPC2EV );
	
	//this->canvas()->resize( this->canvas()->height(), this->canvas()->height() );
	plotLayout()->setAlignCanvasToScales(true);
	emit updated( QString( "PCA Plot" ), 1 );
	replot();
}

void ProRataPCAGraph::MS1selected( long iScan )
{
	if( peppRatioCurrent == NULL )
		return;
	vector< double > vdDataX = peppRatioCurrent->getReferencePeak();
	vector< double > vdDataY = peppRatioCurrent->getTreatmentPeak();
	vector< unsigned long int > viScan = peppRatioCurrent->getScanNumberPeak();
	double dX = 0;
	double dY = 0;
	bool bMS1WithinPeak = false;
	for( unsigned int i = 0; i < viScan.size(); ++i )
	{
	//	qDebug( " in picked %i", (int)viScan[i] );
		if( ( (long)viScan[i] ) == iScan )
		{
			dX = vdDataX[i];
			dY = vdDataY[i];
			bMS1WithinPeak = true;
		}
	}



	if( !bMS1WithinPeak )
		return;

	if( bIsMS1Selected )
		delete qwtpcSelectedMS1;

	qwtpcSelectedMS1 = new QwtPlotCurve( "Selected full scan" );
	bIsMS1Selected = true;
	
	// Create a symbol to represent data  points.
	qwtsMS1Sym.setStyle(QwtSymbol::Ellipse);
	qwtsMS1Sym.setPen(QColor(Qt::green));
	qwtsMS1Sym.setBrush(QBrush(Qt::green));
	qwtsMS1Sym.setSize(6);

	qwtpcSelectedMS1->setSymbol( qwtsMS1Sym );
	qwtpcSelectedMS1->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcSelectedMS1->setYAxis(QwtPlot::yLeft);
	// This curve will not have any connection between points.
	qwtpcSelectedMS1->setStyle( QwtPlotCurve::NoCurve );

	// Add the data to the curve.
	qwtpcSelectedMS1->setData( &dX, &dY, 1 );

	// Add the curve to the graph
	qwtpcSelectedMS1->attach( this );

	// Update the graph.
	replot();
}


