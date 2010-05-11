
#ifndef PRORATAQUANT_H
#define PRORATAQUANT_H

#include <QVariant>
#include <QPushButton>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QFileDialog>
#include <QStringList>

class ProRataQuant : public QWidget
{
	Q_OBJECT

	public:
		ProRataQuant( QWidget* qwParent = 0, Qt::WFlags qwfFl = 0 );
		~ProRataQuant();

		const QStringList & getValues();

	protected:
		void buildUI();
		void setValues();

	protected slots:
		void fastaBrowseSlot();

	private:

		QGroupBox* qgbPeptide;

		QLabel* qtlSmOrd;
		QLabel* qtlLtPkSft;
		QLabel* qtlMinLogRatio;
		QLabel* qtlLogSNRCutff;

		// smoothing order
		QComboBox* qleSmOrd;

		// left peak shift
		QSpinBox* qteLtPkSft;

		// minimum peptide log2ratio
		QSpinBox* qleMinLogRatio;
		
		// log2SNR cutoff
		QDoubleSpinBox* qleLogSNRCutff;

		QLabel* qtlSmWinSize;
		QLabel* qtlRtPkSft;
		QLabel* qtlMaxLogRatio;
		QLabel* qtlRmAmbPeps;

		// smoothing window size
		QComboBox* qleSmWinSize;

		// right peak shift
		QSpinBox* qteRtPkSft;

		// maximum peptide log2ratio
		QSpinBox* qleMaxLogRatio;
		
		// remove ambiguous peptide
		QComboBox* qcbRmAmbPeps;

		QGroupBox* qgbProtein;

		QGroupBox* qgbMinProtein;
		QLabel* qtlMinLogRatioPro;
		QLabel* qtlMinPepNumber;

		// minimum protein log2ratio
		QSpinBox* qleMinLogRatioPro;

		// minimum peptide number
		QSpinBox* qleMinPepNumber;

		QLabel* qtlLogRatioDis;

		// log2ratio discretization
		QComboBox* qleLogRatioDis;

		QGroupBox* qgbMaxProtein;
		QLabel* qtlMaxLogRatioPro;
		QLabel* qtlMaxLogSNR;
		
		// maximum protein log2ratio
		QSpinBox* qleMaxLogRatioPro;
		
		// maximum log2SNR
		QDoubleSpinBox* qleMaxLogSNR;

		QLabel* qtlCIWidth;

		// maximum confidence interval width
		QDoubleSpinBox* qleCIWidth;

		QGroupBox* qgbFasta;
		QLineEdit* qleFasta;
		QPushButton* qpbFasta;
		
		// Layouts

		QVBoxLayout* qvbProRataQuantLayout;

		QHBoxLayout* qgbPeptideLayout;

		QVBoxLayout* qvbPepLabels1;
		QVBoxLayout* qvbPepEntry1;
		QVBoxLayout* qvbPepLabels2;
		QVBoxLayout* qvbPepEntry2;

		QGridLayout* qgbProteinLayout;

		QHBoxLayout* qgbMinProteinLayout;
		QVBoxLayout* qvbLogMinProLabels;
		QVBoxLayout* qvbLogMinProEntrys;

		QHBoxLayout* qhbLgRtDis;


		QVBoxLayout* qgbMaxProteinLayout;
		QHBoxLayout* qhbProMaxLogs;
		QVBoxLayout* qvbLogMaxProLabels;
		QVBoxLayout* qvbLogMaxProEntrys;
		QHBoxLayout* qhbCIWidth;

		QHBoxLayout* qgbFastaLayout;

		QStringList qslValues;
};

#endif //PRORATAQUANT_H
