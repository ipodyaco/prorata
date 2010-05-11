


#ifndef PRORATAGRAPHDIALOG_H
#define PRORATAGRAPHDIALOG_H

#include <QWidget>
#include <QHBoxLayout>
#include <QAction>
#include <QEvent>
#include <QResizeEvent>
#include <QMoveEvent>

#include "proRataGraph.h"

#include <iostream>
using namespace std;

class ProRataGraphDialog : public QWidget
{
	Q_OBJECT;

	public:
		ProRataGraphDialog( QWidget * qwParent = 0, Qt::WFlags flags = 0  );
		~ProRataGraphDialog();

		void setGraph( const QString qsName, ProRataGraph * );
		ProRataGraph *getGraph()
		{	return prgGraph;	}

		QString getGraphName()
		{	return qsGraphName; 	}

		void resizeEvent ( QResizeEvent * event );
		void moveEvent ( QMoveEvent * event );

	protected:
		void closeEvent(QCloseEvent *event);

	signals:
		void closeAction( ProRataGraphDialog * );
		void resizeSignal( const QSize & );
		void moveSignal( const QPoint & );

	private:

		QHBoxLayout *qhblChildMainLayout;
		ProRataGraph * prgGraph;
		QString qsGraphName;

};

#endif //PRORATAGRAPHDIALOG_H

