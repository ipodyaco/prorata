#ifndef PRORATAWORKDIR_H
#define PRORATAWORKDIR_H

#include <QPushButton>
#include <QGroupBox>
#include <QLineEdit>
#include <QLabel>
#include <QLayout>
#include <QFileDialog>
#include <QSizePolicy>
#include <QMessageBox>
#include <QComboBox>

#include "proRataPreProcess.h"
#include "proRataDTASelect.h"
#include "proRataConfigBrowser.h"
#include "proRataRaw2MzXMLBrowser.h"

class ProRataPreProcessWizard;
class ProRataDTASelect;
class ProRataConfigBrowser;
class ProRataRaw2MzXMLBrowser;

class ProRataWorkingDirBrowser : public QWidget
{
	Q_OBJECT

	public:
		ProRataWorkingDirBrowser( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataWorkingDirBrowser();

		const QString & getWorkingDirectory()
		{	return qsWorkingDirectory;	}
		bool isValid();

	public slots:
		void getDirectory();
		void validate();

	protected:
		void buildUI();

		QWidget *qwParent;

		QGroupBox* qgbInputBox;
		QLabel* qlDirectoryLabel;
		QLineEdit* qleDirectoryEntry;
		QPushButton* qpbDirectoryBrowser;

		QHBoxLayout* qhbMainLayout;
		QGridLayout* qgbInputBoxLayout;

		QString qsWorkingDirectory;
		bool bValidity;

		friend class ProRataPreProcessWizard;
		friend class ProRataDTASelect;
		friend class ProRataConfigBrowser;
		friend class ProRataRaw2MzXMLBrowser;

};

#endif // PRORATAWORKDIR_H
