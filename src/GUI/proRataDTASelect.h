#ifndef PRORATADTASELECT_H
#define PRORATADTASELECT_H

#include <QWidget>
#include <QPushButton>
#include <QGroupBox>
#include <QLineEdit>
#include <QLabel>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QFileDialog>
#include <QSizePolicy>


#include "proRataPreProcess.h"
class ProRataPreProcessWizard;

class ProRataDTASelect : public QWidget
{
	Q_OBJECT

	public:
		ProRataDTASelect( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataDTASelect();

		const QString & getResultFile()
		{	return qsResults; 	}
		bool isValid()
		{	return bValidity;	}

	public slots:
		void getDTAResultsFile();
		void validate();

	protected:
		void buildUI();

		QWidget *qwParent;

		QGroupBox* qgbInputBox;
		QLabel* qlResultFile;
		QLineEdit* qleResultFile;
		QPushButton* qpbResultFileBrowser;

		QVBoxLayout* qvbMainLayout;
		QGridLayout* qgbInputBoxLayout;

		QString qsResults;
		bool bValidity;
		friend class ProRataPreProcessWizard;

};

#endif // PRORATADTASELECT_H
