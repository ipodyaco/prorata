#include <QtGui>

#include "proRataPreProcess.h"
#include <iostream>
using namespace std;

ProRataPreProcessWizard::ProRataPreProcessWizard(QWidget *qwParent)
: ProRataPreProcessAbstract(qwParent)
{
	setNumPages(4);
}

QWidget *ProRataPreProcessWizard::createPage(int iIndex)
{
	switch (iIndex) {
		case 0:
			prwdbFirstPage = new ProRataWorkingDirBrowser(this);
			return prwdbFirstPage;
		case 1:
			prdtasSecondPage = new ProRataDTASelect(this);
			return prdtasSecondPage;
		case 2:
			prcbThirdPage = new ProRataConfigBrowser(this);
			return prcbThirdPage;
		case 3:
			prr2mFourthPage = new ProRataRaw2MzXMLBrowser(this);
			return prr2mFourthPage;
	}
	return 0;
}

void ProRataPreProcessWizard::accept()
{

	cout << "Working Directory = " << prwdbFirstPage->getWorkingDirectory().toAscii().data() << endl;
	cout << "DTA results file = " << prdtasSecondPage->getResultFile().toAscii().data() << endl;
	cout << "Configuration file = " << prcbThirdPage->getConfigurationFile().toAscii().data() << endl;
	cout << "Raw data Directory = " << prr2mFourthPage->getRawDataDirectory().toAscii().data() << endl;
	cout << "Conversion Program name = " << prr2mFourthPage->getConversionProgram().toAscii().data() << endl;

	// Do further processing.
	//
	ProRataPreProcessAbstract::accept();
}
