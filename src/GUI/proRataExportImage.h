#include "qwt_plot.h"
#include "qwt_plot_canvas.h"


class ProRataExportImage : public QwtPlotCanvas
{
	public:
		
		ProRataExportImage( QwtPlot *plot );
		~ProRataExportImage(){}

		QPixmap * getPixmap();
};
