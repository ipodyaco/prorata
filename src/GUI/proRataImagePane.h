
#ifndef PRORATAIMAGEPANE_H
#define PRORATAIMAGEPANE_H

#include <QSplitter>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>

#include "proRataImage.h"

#include <vector>

using namespace std;

class ProRataImagePane : public QWidget
{
	Q_OBJECT;

	public:
		ProRataImagePane( QWidget * qwParent = 0  );
		~ProRataImagePane();

		void addImage( ProRataImage * );

	private:
    		QHBoxLayout *qvbMainLayout;
		QSplitter *qsImageDivider;
		vector<ProRataImage *> vpriGraphs;

};

#endif //PRORATAIMAGEPANE_H
