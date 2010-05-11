#ifndef PRORATAPREPROCESS_H
#define PRORATAPREPROCESS_H

#include "proRataPreProcessAbstract.h"
#include "proRataConfigBrowser.h"
#include "proRataDTASelect.h"
#include "proRataRaw2MzXMLBrowser.h"
#include "proRataWorkingDirBrowser.h"

class ProRataWorkingDirBrowser;
class ProRataDTASelect;
class ProRataConfigBrowser;
class ProRataRaw2MzXMLBrowser;

class ProRataPreProcessWizard : public ProRataPreProcessAbstract
{
	Q_OBJECT

	public:
		ProRataPreProcessWizard(QWidget *qwParent = 0);

	protected:
		QWidget *createPage(int iIndex);
		void accept();

	private:
		ProRataWorkingDirBrowser *prwdbFirstPage;
		ProRataDTASelect  *prdtasSecondPage;
		ProRataConfigBrowser  *prcbThirdPage;
		ProRataRaw2MzXMLBrowser  *prr2mFourthPage;

		friend class ProRataWorkingDirBrowser;
		friend class ProRataDTASelect;
		friend class ProRataConfigBrowser;
		friend class ProRataRaw2MzXMLBrowser;

};

#endif // PRORATAPREPROCESS_H
