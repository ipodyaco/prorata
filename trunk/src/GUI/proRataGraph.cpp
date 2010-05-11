
#include "proRataGraph.h"
#include <QPainter>
#include <QMessageBox>
#include "qwt_plot_layout.h"

ProRataGraph::ProRataGraph( QWidget * qwPane ) : QwtPlot( qwPane )
{
	// Create a legend and set proper properties of this legend
	// add it to the graph.
//	qwtlLegend = new QwtLegend;
//	qwtlLegend->setFrameStyle( QFrame::Box|QFrame::Sunken );
//	qwtlLegend->setBackgroundRole( QPalette::Base );
	
	// Future upgrade
	// Make this clickable for showing/hiding the curve
	//qwtlLegend->setItemMode( QwtLegend::ClickableItem );
	
//	insertLegend( qwtlLegend, QwtPlot::BottomLegend );

	// Set proper properties of the canvas.
	setBackgroundRole(QPalette::Light);
	setCanvasBackground( QColor( Qt::white ) );
	canvas()->setFrameStyle(QFrame::Box | QFrame::Plain );

	// Some constant properties.
	canvas()->setLineWidth( 1 );
	setMargin( 1 );
	resize( 400, 400 );

	// Setup the picker
	qwtppPicker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft,
			QwtPicker::PointSelection | QwtPicker::DragSelection, 
			QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn, 
			canvas());
	qwtppPicker->setRubberBandPen( QColor( Qt::darkGray ) );
	qwtppPicker->setTrackerPen( QColor( Qt::darkGray ) );
}

ProRataGraph::~ProRataGraph()
{

}
/*
void ProRataGraph::legendStatus( bool bStatus )
{
	if ( bStatus )
	{
		qwtlLegend->show();
	}
	else
	{
		qwtlLegend->hide();
	}

}
*/
void ProRataGraph::proteinUpdated( ProteinRatio * )
{
	QwtPlot::clear();
	replot();
}

void ProRataGraph::peptideUpdated( PeptideRatio * )
{
	QwtPlot::clear();
	replot();
}

void ProRataGraph::cleanGraph()
{
	QwtPlot::clear();
	replot();
}

void ProRataGraph::saveToFile( const QString &qsFilename, const char * cpFormat )
{
	
     QPixmap *qPPix;
	 QPixmap qpAxesPix = QPixmap::grabWidget( this );
     ProRataExportImage *preiImageCanvas;
     
     preiImageCanvas = (ProRataExportImage *)this->canvas();
     qPPix = preiImageCanvas->getPixmap();

	 QPainter qptPainter( &qpAxesPix );

	 qptPainter.drawPixmap( plotLayout()->canvasRect(), *qPPix, qPPix->rect() );
	 qpAxesPix.save( qsFilename, cpFormat );
}

