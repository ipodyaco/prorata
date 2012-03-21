#ifndef PRORATARAW2MZXML_H
#define PRORATARAW2MZXML_H

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
#include <QMessageBox>
#include <QComboBox>
#include <QDialog>
#include <QFrame>
#include <QProgressBar>
#include <QTimer>
#include <QStringList>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QCoreApplication>


class ProRataRaw2MzXMLBrowser : public QDialog
{
	Q_OBJECT

	public:
		ProRataRaw2MzXMLBrowser( QWidget* parent = 0 );
		ProRataRaw2MzXMLBrowser( QWidget* parent, const QString & qsDir );
		~ProRataRaw2MzXMLBrowser();

		const QString & getRawDataDirectory()
		{	return qsRawDataDirectory;	}
		const QString & getConversionProgram()
		{	return qsConvProgram;	}
		bool isValid();


	public slots:
		void getRawDirectory();
		void getOutputDirectory();
		void validate();
		void validateOuput();
		void convert();
		
	signals:
		void savedDir( const QString & );		

	protected:
		void buildUI();

		QWidget *qwParent;

		QVBoxLayout* qvbClassLayout;

		QLabel* qlHeading;
		QFrame *qfLine0;

		QGroupBox* qgbConvProg;
		QComboBox* qcbConvProg;

		QGroupBox* qgbInputBox;
		QLabel* qlDirectoryLabel;
		QLineEdit* qleDirectoryEntry;
		QPushButton* qpbDirectoryBrowser;

		QGroupBox* qgbOutputBox;
		QLabel* qlOutputDirLabel;
		QLineEdit* qleOutputDirEntry;
		QPushButton* qpbOutputDirectoryBrowser;

		QHBoxLayout* qhbMainLayout;
		QGridLayout* qgbInputBoxLayout;
		QGridLayout* qgbOutputBoxLayout;
		QGridLayout* qgbConvProgLayout;

		QPushButton *qpbCancelButton;
		QPushButton *qpbOkButton;
		QHBoxLayout *qhblButtonLayout;
		QFrame *qfLine1;
	//	QProgressBar *qpgbStatus;


		QString qsRawDataDirectory;
		QString qsOutputDataDirectory;
		QString qsConvProgram;
		bool bValidity;

};

#endif // PRORATARAW2MZXML_H
