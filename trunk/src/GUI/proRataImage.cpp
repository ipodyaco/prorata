
#include "proRataImage.h"

ProRataImage::ProRataImage( const QString & qsFilename, QWidget * qwParent ) :
	QWidget( qwParent )
{
	if ( qsFilename != QString( "" ) )
	{
		qpImage.load(  qsFilename );
	}
	else
	{
		qpImage = QPixmap( notAvailableImage );
	}

	setup();
}


ProRataImage::ProRataImage( const QPixmap & qpImg, QWidget * qwParent ) :
	QWidget( qwParent )
{
	if ( qpImg.isNull() )
	{
		qpImage = QPixmap( notAvailableImage );
	}
	else
	{
		qpImage = qpImg;
	}

	setup();
}

ProRataImage::~ProRataImage()
{
}

void ProRataImage::setup()
{
	// Image Title
	tlImageTitle = new TitleLabel( QString( "Proteomics" ) );
	tlImageTitle->imageHeader();

	// Graph ( Which is a label )
	qlGraph = new QLabel;
	qlGraph->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	qlGraph->setAlignment(Qt::AlignCenter);
	qlGraph->setMinimumSize(300, 300);

	// Image
	qpImage.scaled(qlGraph->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation );

	// Add the image to the Graph label.
	qlGraph->setPixmap( qpImage );


	// Main layout for the widget.
	qvbMainLayout = new QVBoxLayout;
	qvbMainLayout->addWidget( tlImageTitle );
	qvbMainLayout->addWidget( qlGraph );

	/*

	// Help Button
	//QPixmap helpIcon = QPixmap( helpIcon16 );
	//qpbHelp = new QPushButton( QIcon( QPixmap( helpIcon16 ) ), QString( "" ) );
	//qpbHelp = new QPushButton( QIcon( helpIcon ), QString( "" ) );
	//qpbHelp = new QToolButton;
	//qpbHelp->setIcon( helpIcon );
	//qpbHelp->setSize( helpIcon.size() );
	//qpbHelp->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
	//qpbHelp->setFlat( true );
	//qpbHelp->setSize( QIcon( QPixmap( helpIcon32 ) ).size() );
	//qpbHelp->setAlignment( Qt::AlignLeft );

	//connect( qpbHelp, SIGNAL( clicked() ), this, SLOT( helpClicked() ) );

	// New layout to hold spacer and help button.
	//QHBoxLayout *qhbButtomLayout = new QHBoxLayout;
	//qhbButtomLayout->addItem( new QSpacerItem( 50, 10 ) );
	//qhbButtomLayout->addWidget( qpbHelp );

	qvbMainLayout->addItem( new QSpacerItem( 1, 1 ) );
	qvbMainLayout->addWidget( qpbHelp );
	//qvbMainLayout->addLayout( qhbButtomLayout );
	
	*/

	qvbMainLayout->setSpacing( 2 );
	qvbMainLayout->setMargin( 2 );

	this->setLayout( qvbMainLayout );
}

void ProRataImage::setTitle( const QString & qsTitle )
{
	tlImageTitle->setText( qsTitle );
}

void ProRataImage::resizeEvent( QResizeEvent * /*qreNewSize */ )
{

	QSize qszCurrentSize = qpImage.size();

	qszCurrentSize.scale( qlGraph->size(), Qt::KeepAspectRatio);

	if ( ! qlGraph->pixmap()
			|| qszCurrentSize != qlGraph->pixmap()->size() )
	{
		qlGraph->setPixmap( qpImage.scaled(qlGraph->size(), 
					Qt::KeepAspectRatio, Qt::SmoothTransformation ) );
	}

}

