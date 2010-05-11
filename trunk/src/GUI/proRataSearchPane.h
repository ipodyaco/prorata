#ifndef PRORATASEARCHPANE_H
#define PRORATASEARCHPANE_H

#include <QVariant>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCheckBox>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QVariant>
#include <QWidget>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QLabel;
class QLineEdit;
class QPushButton;
class QCheckBox;

class ProRataSearchPane : public QWidget
{
	Q_OBJECT

	public:
		ProRataSearchPane( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataSearchPane();

	signals:
		void searchStringSignal( QString );

	public slots:
		void searchClicked();
		void showAllClicked();

	protected:

		void buildUI();

		QHBoxLayout* findLayout;

		QLineEdit* qleSearchString;
		QPushButton* qpbSearchButton;
		QPushButton* qpbShowAll;

};

#endif // PRORATASEARCHPANE_H

