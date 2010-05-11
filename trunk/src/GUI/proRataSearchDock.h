#ifndef PRORATASEARCHDOCK_H
#define PRORATASEARCHDOCK_H

#include <QVariant>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QCheckBox>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QVariant>
#include <QDockWidget>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QLabel;
class QLineEdit;
class QPushButton;
class QCheckBox;

class ProRataSearchDock : public QDockWidget
{
	Q_OBJECT

	public:
		ProRataSearchDock( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataSearchDock();


	protected:

		void buildUI();

		QHBoxLayout* findLayout;

		QLabel* qlFindLabel;
		QLineEdit* qleSearchString;
		QPushButton* qpbSearchButton;
		QPushButton* qpbShowAll;

};

#endif // PRORATASEARCHDOCK_H

