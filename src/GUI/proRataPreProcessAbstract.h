#ifndef PRORATAPREPROCESSABSTRACT_H
#define PRORATAPREPROCESSABSTRACT_H

#include <QDialog>
#include <QList>
#include <QLayout>
#include <QPushButton>
#include <QLabel>
#include <QFrame>

class ProRataPreProcessAbstract : public QDialog
{
	Q_OBJECT;

	public:
		ProRataPreProcessAbstract(QWidget *qwParent = 0);
		void setButtonEnabled( bool bEnable );

	protected:
		virtual QWidget *createPage( int iIndex ) = 0;
		void setNumPages( int iPages );

	private slots:
		void backButtonClicked();
		void nextButtonClicked();

	protected:
		void switchPage( QWidget *qwOldPage );

		QList<QWidget *> qlstPageHistory;
		int numPages;

		QLabel *qlHeading;
		QFrame *qfLine1;

		QLabel *qlDescription;
		QFrame *qfLine2;

		QFrame *qfLine3;

		QPushButton *qpbHelpButton;
		QPushButton *qpbCancelButton;
		QPushButton *qpbBackButton;
		QPushButton *qpbNextButton;
		QPushButton *qpbFinishButton;

		QHBoxLayout *qhblButtonLayout;
		QVBoxLayout *qvblMainLayout;
};

#endif // PRORATAPREPROCESSABSTRACT_H
