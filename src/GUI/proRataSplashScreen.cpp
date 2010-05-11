
#include "proRataSplashScreen.h"

ProRataSplashScreen::ProRataSplashScreen()
{
	qpImage =  new QPixmap(proRataSplashScreen_xpm);
    proRataSplash = new QSplashScreen(*qpImage);
    proRataSplash->show();
	iTimeInterval = 6;
	
	time(&ttStartTime);
	do
	{
		 time(&ttCurrentTime);
	}
	while((ttCurrentTime - ttStartTime) < iTimeInterval);
}

ProRataSplashScreen::~ProRataSplashScreen()
{

}

void ProRataSplashScreen::execute( QWidget * qwWidget )
{

	proRataSplash->finish( qwWidget );
}
