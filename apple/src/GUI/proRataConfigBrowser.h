#ifndef PRORATACONFIGBROWSER_H
#define PRORATACONFIGBROWSER_H

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

#include "proRataConfigDialog.h"

#include "proRataPreProcess.h"
class ProRataPreProcessWizard;

class ProRataConfigBrowser : public QWidget
{
	Q_OBJECT

	public:
		ProRataConfigBrowser( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataConfigBrowser();

		const QString & getConfigurationFile()
		{	return qsConfigFile;	}
		bool isValid();

	public slots:
		void getConfigFile();
		void newConfigFile();
		void validate();

	protected:
		void buildUI();

		QWidget *qwParent;

		QGroupBox* qgbInputBox;
		QLabel* qlFileLabel;
		QLineEdit* qleFileEntry;
		QPushButton* qpbConfigFileBrowser;
		QPushButton* qpbNew;

		QVBoxLayout* qvbMainLayout;
		QGridLayout* qgbInputBoxLayout;

		QString qsConfigFile;
		bool bValidity;

		friend class ProRataPreProcessWizard;


};

#endif // PRORATACONFIGBROWSER_H
