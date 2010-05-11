
#ifndef PRORATATEXTAREA_H
#define PRORATATEXTAREA_H

#include <QTextBrowser>
#include <QTextCursor>
#include <QFile>
#include <QIODevice>
#include <QString>
#include <string>

#include "proteinInfo.h"
#include "peptideInfo.h"
#include "proteinRatio.h"
#include "peptideRatio.h"

class ProRataTextArea : public QTextBrowser
{
	Q_OBJECT;
	public:
		ProRataTextArea( QWidget * qwParent = 0 );
		~ProRataTextArea();

		void loadContent();


	public slots:

		void substituteContent( QString &, const QString &, const QString & );

		void updateProteinDesc( ProteinInfo * );
		void updatePeptideDesc( PeptideInfo * );
		void updateProteinDesc( ProteinRatio * );
		void updatePeptideDesc( PeptideRatio * );
		void cleanUpPeptide();
		void cleanUp();

	protected:

		QString qsHTMLpepCurrent;
		QString qsHTMLproCurrent;

		QString qsHTMLpepOrigninal;
		QString qsHTMLproOrigninal;

		QString qsHTMLContentCurrent;
		QString qsHTMLContentOriginal;

};
#endif //PRORATATEXTAREA_H
