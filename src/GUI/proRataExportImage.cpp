#include "proRataExportImage.h"

ProRataExportImage::ProRataExportImage( QwtPlot *plot ) : QwtPlotCanvas( plot ) 
{
     
}

QPixmap * ProRataExportImage::getPixmap()
{
     return this->cache();
}
