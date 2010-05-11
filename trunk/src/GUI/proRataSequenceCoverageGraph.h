
#ifndef PRORATASEQUENCECOVERAGEGRAPH_H
#define PRORATASEQUENCECOVERAGEGRAPH_H

#include "proRataGraph.h"
#include "proteinInfo.h"
#include "peptideInfo.h"
#include "qwt_data.h"
#include "fastaData.h"
#include <QString>
#include <QLabel>

class ProRataSequenceCoverageGraph : public ProRataGraph
{
	Q_OBJECT;
	public:
		ProRataSequenceCoverageGraph( QWidget * qwParent = 0  );
		~ProRataSequenceCoverageGraph();

	public slots:
		virtual void proteinUpdated( ProteinInfo * );		
		virtual void cleanGraph();
		virtual void peptideUpdated( PeptideInfo * pPepInfo );
		
	protected:
		FASTAdata * pFASTAdata;
		string sProteinSequence;

		void addPeptide( double dX, int iY1, int iY2 );
		void setXAxisRange( double dLow, double dHigh );
		void setYAxisRange( int iLow, int iHigh );

		QwtPlotCurve *qwtpcSelectedPeptide;
		bool bPreviousSelected;
};
#endif //PRORATASEQUENCECOVERAGEGRAPH_H
