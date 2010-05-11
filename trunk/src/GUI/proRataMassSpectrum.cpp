
#include "proRataMassSpectrum.h"
#include <QMessageBox>
#include <QtGlobal>

#include <iostream>
#include <QVector>
using namespace std;

ProRataMassSpectrum::ProRataMassSpectrum( QWidget * qwPane) : ProRataGraph( qwPane )
{
	// Set the title of the Graph.
	setTitle( "MS1 Scan" );

	// Create a legend and set proper properties of this legend
	// add it to the graph.
//	qwtlLegend = new QwtLegend;
//	qwtlLegend->setFrameStyle( QFrame::Box|QFrame::Sunken );
//	qwtlLegend->setBackgroundRole( QPalette::Base );
//	insertLegend( qwtlLegend, QwtPlot::BottomLegend );
	

	// Set X axis label.
	setAxisTitle(QwtPlot::xBottom, "m/z");
	setAxisTitle(QwtPlot::yLeft, "Relative intensity");
	setAxisLabelFormat( QwtPlot::yLeft, 'f', 0, 3 );

	QwtPlot::replot();
    	qwtZoomer = new QwtPlotZoomer(QwtPlot::xBottom, QwtPlot::yLeft,
			QwtPicker::DragSelection, QwtPicker::AlwaysOff, canvas());
    	qwtZoomer->setRubberBand(QwtPicker::NoRubberBand);
	
	// set the select flag for plot picker
	qwtppPicker->setSelectionFlags(QwtPicker::RectSelection | QwtPicker::DragSelection);
	qwtppPicker->setRubberBand(QwtPicker::RectRubberBand);
  
	// ?? no need to connect picker and zoomer?
//	connect( qwtppPicker, SIGNAL( selected( const QwtDoubleRect & ) ),
	//	qwtZoomer, SLOT( zoom( const QwtDoubleRect & ) ) );

	qwtdFullScanData = NULL;
	ppepActiveRatio = NULL;

	bZoomBaseSet = false;
}

ProRataMassSpectrum::~ProRataMassSpectrum()
{
}

void ProRataMassSpectrum::cleanGraph()
{
	emit updated( QString( "MS1 Scan" ), 0 );
	ProRataGraph::cleanGraph();
}

void ProRataMassSpectrum::setData( ProRataMassSpectrumData * pMassSpecDataInput  )
{
	pMassSpecData = pMassSpecDataInput;
}

void ProRataMassSpectrum::setFullScanData( const QwtArray<double> & dXValues, const QwtArray<double> & dYValues, const QString & qsName )
{
	// Create the data interface
	qwtdFullScanData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
	qwtpcFullScanCurve = new QwtPlotCurve( qsName );
	qwtpcFullScanCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcFullScanCurve->setPen(QPen(Qt::black));
	qwtpcFullScanCurve->setStyle( QwtPlotCurve::Sticks );
	//qwtpcFullScanCurve->setSymbol( qwtsFullScanSym );

	// Set the data to the curve.
	qwtpcFullScanCurve->setData( (*qwtdFullScanData ) );

	// Add the curve to the graph
	qwtpcFullScanCurve->attach( this );

	// Update the graph.
//	replot();
}

void ProRataMassSpectrum::setMZRange( vector< float > vfLower, vector< float > vfUpper, const QColor & qcolor, const QString & qsMZRangeName  )
{
	double xTemp, yTemp;
	if ( ! qwtdFullScanData )
	{
		return;
	}

	int iDataSize = qwtdFullScanData->size();
	int iRangeCount = vfLower.size();

	QwtArray<double> qwtdXData;
	QwtArray<double> qwtdYData;

	int k;
	int n;

	for ( k = 0 ; k < iDataSize; ++k )
	{
		xTemp = qwtdFullScanData->x( k );
		yTemp = qwtdFullScanData->y( k );
		
		for( n = 0; n < iRangeCount; ++n )
		{
			if ( (float)xTemp >= vfLower[n] && (float)xTemp <= vfUpper[n] )
			{
				qwtdXData.push_back( xTemp );
				qwtdYData.push_back( yTemp );
				break;
			}
		}
	}

	QwtArrayData *qwtdNewData = new QwtArrayData( qwtdXData, qwtdYData );

	// Create a curve 
	QwtPlotCurve *qwtNewCurve = new QwtPlotCurve( qsMZRangeName );;
	qwtNewCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtNewCurve->setPen(QPen( qcolor ));
	qwtNewCurve->setStyle( QwtPlotCurve::Sticks );

	// Set the data to the curve.
	qwtNewCurve->setData( (*qwtdNewData ) );

	// Add the curve to the graph
	qwtNewCurve->attach( this );
}


void ProRataMassSpectrum::peptideUpdated( PeptideRatio * prTemp )
{

	if ( !( prTemp ) )
	{
		QwtPlot::clear();
		replot();
		return;
	}

	ppepActiveRatio = prTemp;
	QwtPlot::clear();
	setAxisAutoScale( QwtPlot::xBottom );
	setAxisAutoScale( QwtPlot::yLeft );
	bZoomBaseSet = false;
	setTitle( "MS1 Scan" );
	replot();
}

void ProRataMassSpectrum::updateMSGraph( long iScanInput )
{
//	qDebug( "updataMS %d", iScanInput);
	
	if ( ppepActiveRatio == NULL )
		return;

	if ( pMassSpecData == NULL )
		return;

	unsigned long int iScan = ( unsigned long int )iScanInput;
	
	QwtPlot::clear();

	vector<double> vdMass;
	vector<double> vdRelativeInten;
	
	if( !pMassSpecData->getScan( ppepActiveRatio->getMSfilename(), iScan, vdMass, vdRelativeInten ) )
		return;

	QString qsScanName = QString( "MS1 Scan " ) + QString::number( iScanInput );

	setTitle( qsScanName );
	
	setFullScanData(QVector<double>::fromStdVector( vdMass ),
			QVector<double>::fromStdVector( vdRelativeInten ),  qsScanName  );
	
	vector<float> vfRefLow;
	vector<float> vfRefHi;
	vector<float> vfTretLow;
	vector<float> vfTretHi;

	ppepActiveRatio->getReferenceMZrange( vfRefLow, vfRefHi );
	ppepActiveRatio->getTreatmentMZrange( vfTretLow, vfTretHi );

	// set reference MZ range
	setMZRange( vfRefLow, vfRefHi, QColor( Qt::blue ),
				QString( ProRataConfig::getDenominatorIsotopologue().c_str() )  );

	// set treatment MZ range
	setMZRange( vfTretLow, vfTretHi, QColor( Qt::red ), 
				QString( ProRataConfig::getNumeratorIsotopologue().c_str() ) );
	emit updated( QString( "MS1 Scan" ), 1 );

	replot();

	if( !bZoomBaseSet )
	{
		bZoomBaseSet = true;
		qwtZoomer->setZoomBase();
	}
		
}


