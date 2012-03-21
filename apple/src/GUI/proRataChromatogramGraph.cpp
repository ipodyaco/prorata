
#include "proRataChromatogramGraph.h"
#include <QtGlobal>

ProRataChromatogramGraph::ProRataChromatogramGraph( QWidget * qwPane ) : ProRataGraph( qwPane )
{
	// Set the title of the Graph.
	setTitle( "Selected Ion Chromatograms" );

	qwtlLegend = new QwtLegend;
	qwtlLegend->setFrameStyle( QFrame::Box|QFrame::Sunken );
	qwtlLegend->setBackgroundRole( QPalette::Base );
	insertLegend( qwtlLegend, QwtPlot::BottomLegend );

	bFullScanMarkerExists = false;
	bChroExists = false;

	iMS1ScanNumber = 0;

	qwtpmFullScan = NULL;

	// Set X axis label.
	setAxisTitle(QwtPlot::xBottom, "Retention time (min)");
	setAxisTitle(QwtPlot::yLeft, "Intensity");
	setAxisLabelFormat( QwtPlot::yLeft, 'e', 1, 2 );
	connect( qwtppPicker, SIGNAL( selected( const QwtDoublePoint & ) ),
			this, SLOT( leftPressed( const QwtDoublePoint & ) ) );

	
}

ProRataChromatogramGraph::~ProRataChromatogramGraph()
{
}

void ProRataChromatogramGraph::setStartBoundaryMarker( double dValue )
{

	// Create a marker to show start boundary
	qwtpmStartBoundary = new QwtPlotMarker();
	qwtpmStartBoundary->setLinePen( QPen( QColor( Qt::darkGreen ), 0, Qt::DashDotLine) );
	qwtpmStartBoundary->setLabelAlignment(Qt::AlignLeft|Qt::AlignTop);
	qwtpmStartBoundary->setLineStyle(QwtPlotMarker::VLine);
	//qwtpmStartBoundary->setVisible( false );
	qwtpmStartBoundary->attach( this );
//	qwtpmStartBoundary->setLabel( QString::number( dValue ) );
//	qwtpmStartBoundary->setLabelPen( QPen( QColor( Qt::darkGreen ) ) );
	qwtpmStartBoundary->setXValue( dValue );
//	replot();

}

void ProRataChromatogramGraph::setEndBoundaryMarker( double dValue )
{
	// Create a marker to show end boundary
	qwtpmEndBoundary = new QwtPlotMarker();
	qwtpmEndBoundary->setLinePen( QPen( QColor( Qt::darkGreen ) , 0, Qt::DashDotLine) );
	qwtpmEndBoundary->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmEndBoundary->setLineStyle(QwtPlotMarker::VLine);
	//qwtpmEndBoundary->setVisible( false );
	qwtpmEndBoundary->attach( this );
//	qwtpmEndBoundary->setLabel(  QString::number( dValue ) );
//	qwtpmEndBoundary->setLabelPen( QPen( QColor( Qt::darkGreen ) ) );
	qwtpmEndBoundary->setXValue( dValue );
//	replot();
}

void ProRataChromatogramGraph::setMSMSScanMarker( double dValue )
{
	QwtPlotMarker *qwtpmPPCChro = new QwtPlotMarker;
	// QColor( 192, 131, 0 ) is brown
	qwtpmPPCChro->setLinePen( QPen( QColor( 192, 131, 0) , 0, Qt::DashDotLine) );
	qwtpmPPCChro->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmPPCChro->setLineStyle(QwtPlotMarker::VLine);
//	qwtpmPPCChro->setLabel( QString::number( dValue ) );
//	qwtpmPPCChro->setLabelPen( QPen( QColor( 39, 249, 7) ) );
	qwtpmPPCChro->setXValue( dValue );
	qwtpmPPCChro->attach( this );
//	replot();
}

void ProRataChromatogramGraph::setReferenceData( const QwtArray<double> & dXValues, 
		const QwtArray<double> &dYValues )
{

	// Create the data interface
	qwtdReferenceData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
//	qwtpcReferenceCurve = new QwtPlotCurve( "Reference" );
	qwtpcReferenceCurve = new QwtPlotCurve( QString( ProRataConfig::getDenominatorIsotopologue().c_str() ) );
	qwtpcReferenceCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcReferenceCurve->setPen(QPen(Qt::blue));

	// Set the data to the curve.
	qwtpcReferenceCurve->setData( (*qwtdReferenceData ) );

	// Add the curve to the graph
	qwtpcReferenceCurve->attach( this );

	// Update the graph.
//	replot();
}

void ProRataChromatogramGraph::setTreatmentData( const QwtArray<double> & dXValues, 
		const QwtArray<double> &dYValues )
{

	// Create the data interface
	qwtdTreatmentData = new QwtArrayData( dXValues, dYValues );

	// Create a curve to represent the given data.
//	qwtpcTreatmentCurve = new QwtPlotCurve( "Treatment" );
	qwtpcTreatmentCurve = new QwtPlotCurve( QString( ProRataConfig::getNumeratorIsotopologue().c_str() ) );
	qwtpcTreatmentCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtpcTreatmentCurve->setPen(QPen(Qt::red));

	// Set the data to the curve.
	qwtpcTreatmentCurve->setData( (*qwtdTreatmentData ) );

	// Add the curve to the graph
	qwtpcTreatmentCurve->attach( this );

	// Update the graph.
	replot();
}

void ProRataChromatogramGraph::proteinUpdated( ProteinRatio * )
{
	QwtPlot::clear();
	bFullScanMarkerExists = false;
	bChroExists = false;
	iMS1ScanNumber = 0;

	emit updated( QString( "Ion Chromatogram" ), 0 );

	replot();
}

void ProRataChromatogramGraph::peptideUpdated( PeptideRatio *peppRatio )
{

	if ( !( peppRatio ) )
	{
		QwtPlot::clear();
		bChroExists = false;
		replot();
		return;
	}

	QwtPlot::clear();
	bFullScanMarkerExists = false;
	bChroExists = true;
	iMS1ScanNumber = 0;
	peppRatioCurrent = peppRatio;

	vector<double> vdTemp;
	vector<float> vfTemp = peppRatio->getRetentionTime();

	for( unsigned int i = 0; i < vfTemp.size(); i++ )
	{
		vdTemp.push_back( (double)(vfTemp.at(i)) );
	}

	setReferenceData( 
		QVector<double>::fromStdVector( vdTemp ),
		QVector<double>::fromStdVector( peppRatio->getReferenceChro() ) );

	setTreatmentData( 
		QVector<double>::fromStdVector( vdTemp ),
		QVector<double>::fromStdVector( peppRatio->getTreatmentChro() ) );

	// if the right and left valley time are the same
	// then the peak picking failed and the boundary markers will not be drawn
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
	emit updated( QString( "Ion Chromatogram" ), 1 );
	replot();
}

void ProRataChromatogramGraph::leftPressed( const QwtDoublePoint &pos )
{
//	qDebug(" left pressed " );
	// if there is no chromatogram, don't do anything
	if( !bChroExists )
		return;
	
	iMS1ScanNumber = peppRatioCurrent->getFullScan4Time( pos.x() );
	
	// if there exists a full scan marker, delete it
	if( bFullScanMarkerExists )
		delete qwtpmFullScan;
	// creat a new full scan marker
	double dMSTime = (double)( peppRatioCurrent->getFullScanTime( iMS1ScanNumber ) );
	qwtpmFullScan =  new QwtPlotMarker();
	qwtpmFullScan->setLinePen( QPen( QColor( Qt::green ) , 0, Qt::SolidLine) );
	qwtpmFullScan->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmFullScan->setLineStyle(QwtPlotMarker::VLine);	
//	qwtpmFullScan->setLabel(  ( QString( "MS @ " ) + QString::number( dMSTime ) ) );
	qwtpmFullScan->setLabel(  QString::number( dMSTime ) );
	qwtpmFullScan->setLabelPen( QPen( QColor( Qt::green ) ) );
	qwtpmFullScan->setXValue( dMSTime );
	qwtpmFullScan->attach( this );
	bFullScanMarkerExists = true;
	replot();
	
	// emit the signal for updating Mass spec	
	emit MS1selected( (long) iMS1ScanNumber );

	vector< unsigned long int > viMS2ScanNumber = peppRatioCurrent->getMS2ScanNumber();
	
	if( viMS2ScanNumber.size() > 0 )
	{
		int iMinimumScanDifference = abs((int)viMS2ScanNumber[0] - (int)iMS1ScanNumber );
		unsigned long int iMS2ScanNumber = viMS2ScanNumber[0]; 	
		for( unsigned int i = 1; i < viMS2ScanNumber.size(); ++i )
		{
			// fix the static cast later
			if( abs((int)viMS2ScanNumber[i] - (int)iMS1ScanNumber )  < iMinimumScanDifference )
			{
				iMinimumScanDifference = abs((int)viMS2ScanNumber[i] - (int)iMS1ScanNumber );
				iMS2ScanNumber = viMS2ScanNumber[i];
			}
		}
		if( iMinimumScanDifference < 30 )
			emit MS2selected( (long) iMS2ScanNumber );
	}
	
	
}

/*
 * use scroll wheel to control MS1 selection
 * most code is the same as leftPressed()
 * but not MS2selected signal is emit
 */
void ProRataChromatogramGraph::wheelEvent(QWheelEvent *event)
{
	
	if( !bChroExists )
		return;
	
//	qDebug(" wheel %d ", event->delta() );
	if( event->delta() > 0 )
	{
		iMS1ScanNumber = peppRatioCurrent->getNextFullScan( iMS1ScanNumber, true );
	}
	else
	{
		iMS1ScanNumber = peppRatioCurrent->getNextFullScan( iMS1ScanNumber, false );
	}
	
	// if there exists a full scan marker, delete it
	if( bFullScanMarkerExists )
		delete qwtpmFullScan;
	// creat a new full scan marker
	double dMSTime = (double)( peppRatioCurrent->getFullScanTime( iMS1ScanNumber ) );
	qwtpmFullScan =  new QwtPlotMarker();
	qwtpmFullScan->setLinePen( QPen( QColor( Qt::green ) , 0, Qt::SolidLine) );
	qwtpmFullScan->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
	qwtpmFullScan->setLineStyle(QwtPlotMarker::VLine);	
//	qwtpmFullScan->setLabel(  ( QString( "MS @ " ) + QString::number( dMSTime ) ) );
	qwtpmFullScan->setLabel( QString::number( dMSTime ) );
	qwtpmFullScan->setLabelPen( QPen( QColor( Qt::green ) ) );
	qwtpmFullScan->setXValue( dMSTime );
	qwtpmFullScan->attach( this );
	bFullScanMarkerExists = true;
	replot();
	
	// emit the signal for updating Mass spec	
	emit MS1selected( (long) iMS1ScanNumber );	

}

void ProRataChromatogramGraph::cleanGraph()
{
	emit updated( QString( "Ion Chromatogram" ), 0 );
	ProRataGraph::cleanGraph();
}
