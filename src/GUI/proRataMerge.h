
#ifndef PRORATAMERGE_H
#define PRORATAMERGE_H

#include <QVariant>
#include <QGroupBox>
#include <QListWidget>
#include <QListWidgetItem>
#include <QPushButton>
#include <QLineEdit>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QFrame>
#include <QFileDialog>
#include <QDialog>
#include <QMessageBox>
#include <QFile>
#include <QDir>
#include <QFileInfo>


class ProRataMerge : public QDialog
{
	Q_OBJECT;

	public:
/*		ProRataMerge( QWidget* qwParent = 0, Qt::WFlags qwfFl = 0 ); */
		ProRataMerge( QWidget* qwParent = 0 );
		~ProRataMerge();

		void buildUI();
		void setValues();

	protected:
		QVBoxLayout* qvbProRataMergeLayout;
		QGridLayout* qbgInputLayout;
		QHBoxLayout* qbgOutputLayout;

		QGroupBox* qbgInput;
		QListWidget* qlwInput;
		QPushButton* qpbAddDir;
		QPushButton* qpbRemoveDir;
		QGroupBox* qbgOutput;
		QLineEdit* qleOutput;
		QPushButton* qpbOuputBrowse;


		QPushButton *qpbCancelButton;
		QPushButton *qpbMergeButton;
		QHBoxLayout *qhblButtonLayout;
		QFrame *qfLine1;

		QString qsOutputDirectory;


	protected slots:
		void merge();
		void addDirectory();
		void removeDirectory();
		void browseForOutputDirectory();
		void validateOutputDir( const QString & qsDir);
};

#endif // PRORATAMERGE_H

