
#ifndef XMLPROCESSOR_H
#define XMLPROCESSOR_H

#include <QtXml>

#include <QDomDocument>
#include <QDomNode>
#include <QDomElement>
#include <QString>
#include <QFile>

#include "proRataProteinTuple.h"
#include "proRataPeptideTuple.h"

#include <iostream>

using namespace std;

class ProRataXmlProcessors
{
	public:
		ProRataXmlProcessors( const QString & qsFilename );
		ProRataXmlProcessors();
		~ProRataXmlProcessors();


		void setFilename( const QString & qsFilename );

		bool processProtein( QDomElement &qdeProtein );
		bool processPeptide( QDomElement &qdePeptide );

		ProRataProteinTuple * getNextProtein();
		ProRataPeptideTuple * getNextPeptide();

		ProRataProteinTuple * getProtein( const QString & qsLocus );
		ProRataPeptideTuple * getPeptide( const QString & qsLocus, 
				const QString & qsChromatogramId );

		QString getText( QDomElement &qdeTag, const QString & qsTagName );

		void setProtein( const QString & qsLocus );
		void setPeptide( const QString & qsLocus, 
				const QString & qsChromatogramId );

		void resetPointers();

	private:
		
		bool processFile( const QString & qsFilename );

		void goToProteinNode( const QString & qsLocus );
		void goToPeptideNode( const QString & qsChromatogramId );

		QDomDocument qddResults;
		QDomElement qdeRoot;
		QDomElement qdeFirstProteinElement;
		QDomElement qdeFirstPeptideElement;
		QDomElement qdeCurrentProteinElement;
		QDomElement qdeCurrentPeptideElement;
		QDomNode qdnSelectedProtein;
		QDomNode qdnSelectedPeptide;

		ProRataProteinTuple * prptProteinSeletedData;
		ProRataPeptideTuple * prdtPeptideSeletedData;

		QString qsFilename;

};

#endif //XMLPROCESSOR_H
