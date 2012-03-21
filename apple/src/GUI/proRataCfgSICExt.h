
#ifndef PRORATACFGSIC_H
#define PRORATACFGSIC_H

#include <QVariant>
#include <QWidget>
#include <QPushButton>
#include <QGroupBox>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QTextEdit>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QFileDialog>

class ProRataCfgSIC : public QWidget
{
	Q_OBJECT

	public:
		ProRataCfgSIC( QWidget* parent = 0, Qt::WFlags fl = 0 );
		~ProRataCfgSIC();

		void buildUI();
		void setValues();
		const QStringList & getValues();

	protected slots:
		//void existingBrowseSlot();
		//void addCompositionSlot();
		//
		//
		//
		//void fastaBrowseSlot();

	protected:

		QString convertToXML( QString content, QString tag );
		QVBoxLayout* qvbMainLayout;
		QVBoxLayout* qgbNewCfgLayout;
		QHBoxLayout* qlyMsIdFastaLayout;
		QHBoxLayout* qgbMSFileLayout;
		QHBoxLayout* qgbIDFileLayout;
		//QHBoxLayout* qgbFastaLayout;
		//
		QHBoxLayout* qgbRetentionIntLayout;
		QVBoxLayout* qvbRetentionLabelsLayout;
		QVBoxLayout* qvbRetentionEntrysLayout;

		QHBoxLayout* qgbMassToChargeLayout;
		QVBoxLayout* qvbMassToChargeLabelsLayout;
		QVBoxLayout* qvbMassToChargeEntrysLayout;

		QVBoxLayout* qgbConpositionLayout;
		QHBoxLayout* qlyEnrichmentLayout;

		QGroupBox* qgbNewCfg;
		QGroupBox* qgbMSFile;
		QComboBox* qcbMSFile;
		QGroupBox* qgbIDFile;
		QComboBox* qcbIDFile;
		//QGroupBox* qgbFasta;
		//QLineEdit* qleFasta;
		//QPushButton* qpbFasta;
		QGroupBox* qgbRetentionInt;
		QLabel* qlBeforeMS;
		QDoubleSpinBox* qleBeforeMS;
		QLabel* qlAfterMS;
		QDoubleSpinBox* qleAfterMS;
		QLabel* qlBetMS;
		QDoubleSpinBox* qleBetMS;
		QGroupBox* qgbMassToCharge;
		QLabel* qlPlus;
		QDoubleSpinBox* qlePlus;
		QLabel* qlMinus;
		QDoubleSpinBox* qleMinus;
		QLabel* qlCutoff;
		
		QDoubleSpinBox* qleCutff;
		QGroupBox* qgbConposition;
		QTextEdit* qteComposition;

		// enrichment percentage
		QLabel* qlEnrichmentPctg;
		QDoubleSpinBox* qleEnrichmentPctg;

		QStringList qslValues;

};

#endif // PRORATACFGSIC_H
