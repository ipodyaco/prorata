#include "proRataTandemMassSpectrum.h"

ProRataTandemMassSpectrum::ProRataTandemMassSpectrum( QWidget * qwPane ) : ProRataMassSpectrum( qwPane )
{
	setTitle( "MS2 Scan" );
	bZoomBaseSet = false;
	// creat isotopologue objects
	residueMap mAtomicComposition;
	ProRataConfig::getResidueAtomicComposition( mAtomicComposition );
	residueMap::const_iterator iterResidueMap;
	for( iterResidueMap = mAtomicComposition.begin(); iterResidueMap != mAtomicComposition.end(); ++iterResidueMap )
	{
		vpIsotopologue.push_back( new Isotopologue( iterResidueMap->first, iterResidueMap->second ) );
	}

}

ProRataTandemMassSpectrum::~ProRataTandemMassSpectrum()
{

}

void ProRataTandemMassSpectrum::cleanGraph()
{
	emit updated( QString( "MS2 Scan" ), 0 );
	ProRataGraph::cleanGraph();
}


void ProRataTandemMassSpectrum::peptideUpdated( PeptideRatio * prTemp )
{

	if ( !( prTemp ) )
	{
		QwtPlot::clear();
		replot();
		return;
	}

	ppepActiveRatio = prTemp;
	viMS2ScanNumber = ppepActiveRatio->getMS2ScanNumber();
	sSequence = ppepActiveRatio->getSequence();
	iChargeState =  ppepActiveRatio->getChargeState();
	QwtPlot::clear();
	setAxisAutoScale( QwtPlot::xBottom );
	setAxisAutoScale( QwtPlot::yLeft );
	bZoomBaseSet = false;
	setTitle( "MS2 Scan" );

	emit updated( QString( "MS2 Scan" ), 0 );


	replot();
}


void ProRataTandemMassSpectrum::updateMSGraph( long iScanInput )
{
	
	if ( ppepActiveRatio == NULL )
		return;

	if ( pMassSpecData == NULL )
		return;

	unsigned long int iScan = ( unsigned long int )iScanInput;

	int iMinDistance = 10000;
	int iDistance = 0;
	unsigned long int iMS2Scan = 0;
	unsigned int i;
	for( i = 0; i < viMS2ScanNumber.size(); ++i )
	{
		iDistance = abs( (int)(viMS2ScanNumber[i] - iScan) );
		if( iMinDistance > iDistance )
		{
			iMinDistance = iDistance;
			iMS2Scan = viMS2ScanNumber[i];
		}
	}
	
	if( iMinDistance > 25 )
		return;
	
	QwtPlot::clear();

	vector<double> vdMass;
	vector<double> vdRelativeInten;
	double dPrecursorMZ;
	
	if( !pMassSpecData->getScan( ppepActiveRatio->getMSfilename(), iMS2Scan, vdMass, vdRelativeInten, dPrecursorMZ ) )
		return;

	QString qsScanName = QString( "MS2 Scan " ) + QString::number( iMS2Scan ) + " (Precursor " + QString::number( dPrecursorMZ ) + ")";

	setTitle( qsScanName );
	
	setFullScanData(QVector<double>::fromStdVector( vdMass ),
			QVector<double>::fromStdVector( vdRelativeInten ),  qsScanName  );

	double dPeptideMZ = 0;
	double dMinMZerror = 10000;
	int iIsoIndex = 0;
	
	for( i = 0; i < vpIsotopologue.size(); ++i )
	{
		dPeptideMZ = vpIsotopologue[i]->computeAverageMass( sSequence );
		if( dPeptideMZ <= 0.05 )
			return;
		dPeptideMZ = ( dPeptideMZ + (double)iChargeState ) / (double)iChargeState;
		if( dMinMZerror > fabs( dPrecursorMZ - dPeptideMZ ) )
		{
			dMinMZerror = fabs( dPrecursorMZ - dPeptideMZ );
			iIsoIndex = i;
		}
	}	

	vector< double > vdYion;
	vector< double > vdBion;

	vpIsotopologue.at(iIsoIndex)->computeProductIonMass( sSequence, vdYion, vdBion );

	markIonSeries( vdYion, QColor( Qt::blue ), QString( "y" ) );
	markIonSeries( vdBion, QColor( Qt::red ), QString( "b" ) );

	emit updated( QString( "MS2 Scan" ), 1 );

	replot();

	if( !bZoomBaseSet )
	{
		bZoomBaseSet = true;
		qwtZoomer->setZoomBase();
	}
		
}

void ProRataTandemMassSpectrum::markIonSeries( vector< double > vdIonMZ, const QColor & qcolor, const QString & qsName )
{

	double xTemp, yTemp;
	if ( ! qwtdFullScanData )
	{
		return;
	}

	unsigned int iDataSize = qwtdFullScanData->size();
	unsigned int k;
	unsigned int n;
	QwtArray<double> qwtdXData;
	QwtArray<double> qwtdYData;
	
	for(  n = 0; n <  vdIonMZ.size(); ++n )
	{
		double dXlabel = vdIonMZ[n];
		double dYlabel = 0;
		double fLow = vdIonMZ[n] - 1.5;
		double fHi = vdIonMZ[n] + 1.5;

		for( k = 0 ; k < iDataSize; ++k )
		{
			xTemp = qwtdFullScanData->x( k );
			yTemp = qwtdFullScanData->y( k );
			if ( (double)xTemp >= fLow && (double)xTemp <= fHi )
			{
				qwtdXData.push_back( xTemp );
				qwtdYData.push_back( yTemp );
				if( yTemp > dYlabel )
					dYlabel = yTemp;
			}

		}

		if( dYlabel < 3 )
			continue;

		if( dYlabel > 98 )
		{
			dYlabel = 97;
			dXlabel = dXlabel + 18 ;
		}
		
		QwtPlotMarker *qwtpmIon = new QwtPlotMarker;
		qwtpmIon->setLabelAlignment(Qt::AlignCenter|Qt::AlignTop);
		qwtpmIon->setLineStyle(QwtPlotMarker::NoLine);
		qwtpmIon->setLabel( qsName + QString::number( (n+1) ) );
		qwtpmIon->setLabelPen( QPen( qcolor ) );
		qwtpmIon->setXValue( dXlabel );
		qwtpmIon->setYValue( dYlabel + 0 );
		qwtpmIon->attach( this );	
	}

	QwtArrayData *qwtdNewData = new QwtArrayData( qwtdXData, qwtdYData );

	// Create a curve 
	QwtPlotCurve *qwtNewCurve = new QwtPlotCurve(qsName + QString(" Ion Series") );
	qwtNewCurve->setRenderHint(QwtPlotItem::RenderAntialiased);
	qwtNewCurve->setPen(QPen( qcolor ));
	qwtNewCurve->setStyle( QwtPlotCurve::Sticks );

	// Set the data to the curve.
	qwtNewCurve->setData( (*qwtdNewData ) );

	// Add the curve to the graph
	qwtNewCurve->attach( this );
	
}

