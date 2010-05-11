
#include "proRataMainWindow.h"
#include "proRataSplashScreen.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

//	ProRataSplashScreen prssScreen;

	ProRataMainWindow main;
	main.show();

//	prssScreen.execute( &main );
	
	return app.exec();
}

