#ifndef PRORATASPLASHSCREEN_H
#define PRORATASPLASHSCREEN_H

#include <QPixmap>
#include <QSplashScreen>

#include <time.h>

#include "proRataSplash.h"

class ProRataSplashScreen
{

	public:
		ProRataSplashScreen();
		~ProRataSplashScreen();
		void setInterval( int iInt = 6 )
		{	iTimeInterval = iInt;	}
		void execute( QWidget *qwWidget );

	protected:
		QPixmap *qpImage;
		QSplashScreen *proRataSplash;
		time_t ttStartTime, ttCurrentTime;
		int iTimeInterval;

};

#endif  //PRORATASPLASHSCREEN_H
