#ifndef PRORATAEXEC_H
#define PRORATAEXEC_H

#include <QWidget>
#include <QVariant>
#include <QPushButton>
#include <QLabel>
#include <QGroupBox>
#include <QLineEdit>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QFileDialog>
#include <QSizePolicy>
#include <QMessageBox>
#include <QFile>
#include <QDir>
#include <QProcess>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QLabel;
class QGroupBox;
class QLineEdit;
class QPushButton;

class ProRataExec : public QWidget
{
	Q_OBJECT

	public:
		ProRataExec( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataExec();

	public slots:
		void workingDirBrowseSlot();
		void idFileBrowseSlot();
		void cancelSlot();
		void okSlot();


	protected:

			void buildUI();
			QVBoxLayout* qvbMainLayout;
			QGridLayout* qgbInputBoxLayout;
			QHBoxLayout* qhblControls;
			QSpacerItem* spacer1;

			QLabel* qlHeading;
			QGroupBox* qgbInputBox;
			QLineEdit* qleWorkingDir;
			QPushButton* qpbWorkingDirBrowse;
			QLabel* qlIdFile;
			QLineEdit* qleIdFile;
			QPushButton* qpbIdFileBrowse;
			QLabel* qlWorkingDir;
			QPushButton* qpbCancel;
			QPushButton* qpbExec;

};

#endif // PRORATAEXEC_H
