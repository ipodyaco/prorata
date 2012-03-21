
#ifndef TITLELABEL_H
#define TITLELABEL_H

#include <QLabel>
#include <QWidget>
#include <QString>

class TitleLabel : public QLabel
{
	public:
		TitleLabel( const QString & text, QWidget * qwParent = 0 );
		~TitleLabel();

		void tableHeader();
		void imageHeader();

		/*
	protected:
		void resizeEvent( QResizeEvent * );
		*/
};
#endif //TITLELABEL_H
