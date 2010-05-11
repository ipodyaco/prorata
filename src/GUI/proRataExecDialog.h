
#ifndef PRORATAEXECDIALOG_H
#define PRORATAEXECDIALOG_H

#include <QGroupBox>
#include <QLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QDialog>

#include <QString>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QCoreApplication>
#include "proRataConfigDialog.h"
#include "proRataRaw2MzXMLBrowser.h"
#include "proRataMerge.h"

#include <unistd.h>

class QAssistantClient;



class ProRataExecDialog : public QDialog
{
	Q_OBJECT;

	public:
		ProRataExecDialog( QWidget* qwParent = 0 );
		~ProRataExecDialog();

		const QString & getWorkingDirectory()
		{	return qsWorkingDirectory;	}
		//bool isValid();

	public slots:

		void helpExec();
		
		void getWdirDirectory();
		void getDTAResultsFile();
		void getConfigFile();
		void getMZXMLDirectory();
		void creatConfigFile();

		void creatMzxml();
		void creatDTASelect();
	//	void proceedCheck();
		void proceed();
		bool validateWDir(const QString & );
		
	signals:
		void qprFilename( QString );

	protected:
		void buildUI();
		bool validateDTAFile();
		bool validateConfigFile();
		bool validateMZMLDir();

		QVBoxLayout* qvbMainLayout;

		QGroupBox* qgbWDirInputBox;
		QGridLayout* qgbWDirInputBoxLayout;
		QLabel* qlDirectoryLabel;
		QLineEdit* qleDirectoryEntry;
		QPushButton* qpbDirectoryBrowser;



		QGroupBox* qgbDTAInputBox;
		QGridLayout* qgbDTAInputBoxLayout;
		QLabel* qlResultFile;
		QLineEdit* qleResultFile;
		QPushButton* qpbResultFileBrowser;
		QPushButton* qpbCreatDTASelect;



		QGroupBox* qgbCFGInputBox;
		QGridLayout* qgbCFGInputBoxLayout;
		QLabel* qlFileLabel;
		QLineEdit* qleFileEntry;
		QPushButton* qpbConfigFileBrowser;
		QPushButton* qpbNew;


		QGroupBox* qgbMZXMLInputBox;
		QGridLayout* qgbMZXMLInputBoxLayout;
		QLabel* qlMZXMLDirectoryLabel;
		QLineEdit* qleMZXMLDirectoryEntry;
		QPushButton* qpbMZXMLDirectoryBrowser;

		QPushButton *qpbHelpButton;
		QLabel *qlStatus;
		QPushButton *qpbCancelButton;
		QPushButton *qpbMzXMLButton;
		QPushButton *qpbProceedButton;
		QHBoxLayout *qhblButtonLayout;
		QLabel *qlHeading;
	//	QFrame *qfLine1;
	//	QLabel *qlDescription;
	//	QFrame *qfLine2;
		QFrame *qfLine3;

		QString qsWorkingDirectory;
		QString qsDTAResultsFile;
		QString qsCFGFile;
		QString qsMZXMLDir;
};

#endif // PRORATAEXECDIALOG_H

