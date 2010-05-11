
#ifndef PRORATARESIDUECOMPOSITION_H
#define PRORATARESIDUECOMPOSITION_H

#include "proRataIsotop.h"

class ProRataResidueComposition : public QWidget
{
	Q_OBJECT

	public:
		
		ProRataResidueComposition( QWidget* qwParent = 0, Qt::WFlags qwfFl = 0 );
		~ProRataResidueComposition();

		void buildUI();
		void setValues();

	protected:

		ProRataIsotopologue *priReference;
		ProRataIsotopologue *priTreatment;

		QGridLayout* qgMainLayout;

};


#endif //PRORATARESIDUECOMPOSITION_H
