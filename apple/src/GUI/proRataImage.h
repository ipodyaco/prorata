
#ifndef PRORATAIMAGE_H
#define PRORATAIMAGE_H

#include <QString>
#include <QPixmap>
#include <QLabel>
#include <QVBoxLayout>
#include <QSize>
#include <QToolButton>
#include <QIcon>

#if 0
#include "helpIcons.h"
#endif

#include "notAvailableImage.h"

#include "titleLabel.h"

class ProRataImage : public QWidget
{
	public:
		ProRataImage( const QString & qsFilename, QWidget * qwParent = 0 );
		ProRataImage( const QPixmap & qpImg, QWidget * qwParent = 0 );
		~ProRataImage();

		void setTitle( const QString & );

	protected:
		void resizeEvent( QResizeEvent * );

	private slots:
		void helpClicked();

	private:
		void setup();

		QPixmap qpImage;
		QLabel *qlGraph;
		TitleLabel *tlImageTitle;
		QVBoxLayout *qvbMainLayout;
		//QToolButton *qpbHelp;

};
#endif //PRORATAIMAGE_H
