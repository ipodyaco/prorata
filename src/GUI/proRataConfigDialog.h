
#ifndef PRORATACONFIGDIALOG_H
#define PRORATACONFIGDIALOG_H

#include "proRataCfgSICExt.h"
#include "proRataIsotop.h"
#include "proRataQuant.h"
#include "proRataConfig.h"

#include <QLabel>
#include <QTabWidget>
#include <QDir>
#include <QFile>
#include <QtXml>
#include <QMessageBox>
#include <string>
#include <vector>

using namespace std;

class ProRataConfigDialog : public QWidget
{
	Q_OBJECT

	public:
		
		ProRataConfigDialog( QWidget* qwParent = 0, Qt::WFlags qwfFl = 0 );
		~ProRataConfigDialog();

		void buildUI();
		void setValues();
		QWidget * getIsotopologue( const QString & qsName );
		void setConfigVersion( const QString & qsVer = QString::fromStdString( ProRataConfig::getProRataVersion() )  )
		{	qsVersion = qsVer;	}

	protected slots:
	//	void save();
		void saveAs();
	
	signals:
		void savedFilename( const QString & );
		
	protected:

		void buildXML( const QString & );

		void addElement( QDomElement & qdeParent, const QString &qsTagName, 
				const QString &qsValue );

		void addElement( QDomElement & qdeParent, const QString &qsTagName );
		void addComment( QDomElement & qdeParent, const QString &qsValue );

		QTabWidget *qtwMainTab;

		ProRataCfgSIC *prCfgSic;
		ProRataQuant *prQuant;

		QHBoxLayout *qhbControlsLayout;

		QPushButton *qpbSave;
		QPushButton *qpbSaveAs;
		QPushButton *qpbCancel;

		QGridLayout* qgMainLayout;

		QDomDocument qddXmlDoc;
		QString qsVersion;

};

#endif //PRORATACONFIGDIALOG_H


