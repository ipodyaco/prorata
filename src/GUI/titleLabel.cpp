
#include "titleLabel.h"

TitleLabel::TitleLabel( const QString & text, QWidget * qwParent ) :
	QLabel( text, qwParent )
{
	setText( text );
	setAlignment( Qt::AlignHCenter );
	setMargin( 4 );
//	setMinimumSize( size() );
	setIndent(2);
}

TitleLabel::~TitleLabel()
{

}

void TitleLabel::tableHeader()
{
	setText( "<center>" + text() + "<center>" );
}

void TitleLabel::imageHeader()
{
	setText( "<center><b>" + text() + "</b><center>" );
}

/*
void TitleLabel::resizeEvent( QResizeEvent * )
{
	setMinimumSize( size() );
}
*/
