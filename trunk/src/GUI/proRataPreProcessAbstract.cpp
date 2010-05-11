#include <QtGui>

#include "proRataPreProcessAbstract.h"

ProRataPreProcessAbstract::ProRataPreProcessAbstract(QWidget *qwParent)
	: QDialog(qwParent)
{

	qlHeading = new QLabel( tr("<center><b>ProRata data pre-processing</b></center>" ) );
	qfLine1 = new QFrame;
	qfLine1->setLineWidth( 1 );
	qfLine1->setFrameStyle( QFrame::HLine | QFrame::Sunken );

	qlDescription = new QLabel( tr("This wizard will guide you through the ProRata "
				"data pre-processing stages. \nPlease provide the necessary "
				"inputs at various stages. For more information,\nplease click "
			       	"help."	) );
	qfLine2 = new QFrame;
	qfLine2->setLineWidth( 1 );
	qfLine2->setFrameStyle( QFrame::HLine | QFrame::Sunken );


	qfLine3 = new QFrame;
	qfLine3->setLineWidth( 1 );
	qfLine3->setFrameStyle( QFrame::HLine | QFrame::Sunken );

	qpbHelpButton = new QPushButton(tr("&Help"));
	qpbCancelButton = new QPushButton(tr("Cancel"));
	qpbBackButton = new QPushButton(tr("< &Back"));
	qpbNextButton = new QPushButton(tr("Next >"));
	qpbNextButton->setEnabled(false);
	qpbFinishButton = new QPushButton(tr("&Finish"));
	qpbFinishButton->setEnabled(false);

	connect(qpbCancelButton, SIGNAL(clicked()), this, SLOT(reject()));
	connect(qpbBackButton, SIGNAL(clicked()), this, SLOT(backButtonClicked()));
	connect(qpbNextButton, SIGNAL(clicked()), this, SLOT(nextButtonClicked()));
	connect(qpbFinishButton, SIGNAL(clicked()), this, SLOT(accept()));

	qhblButtonLayout = new QHBoxLayout;
	qhblButtonLayout->addWidget(qpbHelpButton);
	qhblButtonLayout->addStretch(1);
	qhblButtonLayout->addWidget(qpbCancelButton);
	qhblButtonLayout->addWidget(qpbBackButton);
	qhblButtonLayout->addWidget(qpbNextButton);
	qhblButtonLayout->addWidget(qpbFinishButton);

	qvblMainLayout = new QVBoxLayout;
	qvblMainLayout->addWidget( qlHeading );
	qvblMainLayout->addWidget( qfLine1 );
	qvblMainLayout->addWidget( qlDescription );
	qvblMainLayout->addWidget( qfLine2 );
	qvblMainLayout->addWidget( qfLine3 );
	qvblMainLayout->addLayout(qhblButtonLayout);
	setLayout(qvblMainLayout);
}

void ProRataPreProcessAbstract::setButtonEnabled(bool bEnable)
{
	if (qlstPageHistory.size() == numPages)
		qpbFinishButton->setEnabled(bEnable);
	else
		qpbNextButton->setEnabled(bEnable);
}

void ProRataPreProcessAbstract::setNumPages( int iPages )
{
	numPages = iPages;
	qlstPageHistory.append(createPage(0));
	switchPage(0);
	setFixedWidth( sizeHint().width() );
}

void ProRataPreProcessAbstract::backButtonClicked()
{
	qpbNextButton->setEnabled(true);
	qpbFinishButton->setEnabled(true);

	QWidget *qwOldPage = qlstPageHistory.takeLast();
	switchPage(qwOldPage);
	delete qwOldPage;
}

void ProRataPreProcessAbstract::nextButtonClicked()
{
	//qpbNextButton->setEnabled(true);
	qpbNextButton->setEnabled(false);
	//qpbFinishButton->setEnabled(qlstPageHistory.size() == numPages - 1);

	QWidget *qwOldPage = qlstPageHistory.last();
	qlstPageHistory.append(createPage(qlstPageHistory.size()));
	switchPage(qwOldPage);
}

void ProRataPreProcessAbstract::switchPage(QWidget *qwOldPage)
{
	if (qwOldPage) {
		qwOldPage->hide();
		qvblMainLayout->removeWidget(qwOldPage);
	}

	QWidget *newPage = qlstPageHistory.last();
	qvblMainLayout->insertWidget(4, newPage);
	newPage->show();
	newPage->setFocus();

	qpbBackButton->setEnabled(qlstPageHistory.size() != 1);
	if (qlstPageHistory.size() == numPages) {
		qpbNextButton->setEnabled(false);
		//qpbFinishButton->setDefault(true);
	} 
	else 
	{
		qpbNextButton->setDefault(true);
		qpbFinishButton->setEnabled(false);
	}

	setWindowTitle(tr("ProRata Pre-processing - Step %1 of %2")
			.arg(qlstPageHistory.size())
			.arg(numPages));
}
