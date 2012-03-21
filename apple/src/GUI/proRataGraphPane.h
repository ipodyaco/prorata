
#ifndef PRORATAGRAPHPANE_H
#define PRORATAGRAPHPANE_H

#include <QSplitter>
#include <QWidget>
#include <QVBoxLayout>
#include <QTabWidget>
#include <QMenu>
#include <QAction>
#include <QContextMenuEvent>
#include <QDialog>
#include <QPoint>
#include <QEvent>

#include "proRataGraph.h"
#include "proRataGraphDialog.h"

#include <vector>
#include <iostream>
using namespace std;

class ProRataGraphPane : public QWidget
{
	Q_OBJECT;

	public:
		ProRataGraphPane( QWidget * qwParent = 0  );
		~ProRataGraphPane();

		void addGraph( const QString qsName, ProRataGraph * );

		/*
		void setProteinPane( ProRataImagePane *pane ){
			pripProteinPane = pane;
		}

		void setPeptidePane( ProRataImagePane *pane ){
			pripPeptidePane = pane;
		}
		*/

	protected:
		virtual void contextMenuEvent(QContextMenuEvent *event);

	private slots:
		void detach();
		void exportGraph();
		void closeDialog( ProRataGraphDialog * );
		void resizeSlot( const QSize &qszSz );
		void moveSlot( const QPoint &qptPt );
		void print();
	//	void graphUpdated(QString, int);
	//	void currentSelected(int);

	private:

		void buildMenu();
		//ProRataGraph * getGraph( const QString & );

    		QVBoxLayout *qvbMainLayout;
		vector< ProRataGraph *> vprgGraphs;

		QTabWidget *qtwGraphTab;

		QMenu *qmContextMenu;
		
		/*
		QAction *qaDetachLikelihood;
		QAction *qaDetachChromatogram;
		QAction *qaDetachPPC;
		QAction *qaDetachMassSpectrum;
		QAction *qaDetachPCA;
		QAction *qaDetachSequenceCoverage;
		*/

		QAction *qaDetach;
		QAction *qaExport;
		//QAction *qaPrintAll;

		QSize qszDialogSize;
		QPoint qpDialogPosition;
	//	QIcon qiUpdated;
	//	QIcon qiSeen;
};

#endif
