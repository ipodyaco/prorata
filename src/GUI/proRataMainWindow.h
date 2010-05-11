#ifndef PRORATAMAINWINDOW_H
#define PRORATAMAINWINDOW_H

#include <QMainWindow>
#include <QtGui>
#include <QAction>
#include <QMenu>
#include <QTextEdit>


#include <iostream>

using namespace std;

#include "proRataMainWidget.h"
#include "proRataExec.h"
#include "helpIcons.h"
#include "proRataMainIcon.h"

//#include "proRataSearchDock.h"
//#include "proRataConfigDialog.h"
#include "proRataPreProcess.h"
#include "proRataExecDialog.h"
#include "proRataRaw2MzXMLBrowser.h"
#include "proRataMerge.h"

class ProRataMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    ProRataMainWindow(QWidget *parent = 0, Qt::WFlags flags = 0);

/*
protected:
    void closeEvent(QCloseEvent *event);
*/
private slots:

    //File
    void open();
    void open( QString );
    //void openRecent();
    void print();
    void mainClose();
    void quit();

    /*
    //Edit
    void cut();
    void copy();
    void paste();
    */

    //View
    void defaultView();
    void search();
	void showTables();
    void showText();
    void showGraphs();

    //Tools
    void executeProRata();
    void configuration();
	void mzXMLConversion();
	void mergeDirId();
	
    //void wizard();
    //void generateMZXML();

    //Help
    void contents();
    void generalHelp();
    void introduction();
    void about();


private:
    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();
    void readSettings();

    void loadFile(const QString &qsFileName);
    void setCurrentFile(const QString &qsFileName);
    QString strippedName(const QString &qsFullFileName);

    QMenu *qmFileMenu;
    // QMenu *qmEditMenu;
    QMenu *qmViewMenu;
    QMenu *qmToolsMenu;
    QMenu *qmHelpMenu;
    
    QAction *qaOpenAct;
    //QAction *qaOpenRecentAct;
    QAction *qaPrintAct;
    QAction *qaCloseAct;
    QAction *qaQuitAct;

    /*
    QAction *qaCutAct;
    QAction *qaCopyAct;
    QAction *qaPasteAct;
    */

    QAction *qaDefaultViewAct;
    QAction *qaSearch;
	QAction *qaShowTablesAct;
    QAction *qaShowTextAct;
    QAction *qaShowGraphsAct;

    QAction *qaExecuteProRata;
	QAction *qaConfiguration;
	QAction *qaMzXMLConvert;
	QAction *qaMerge;

    //QAction *qaWizardAct;
    //QAction *qaGenerateMZXMLAct;
    //QAction *qaOptionsAct;

    QAction *qaContentsAct;
    QAction *qaGeneralHelpAct;
    QAction *qaIntroductionAct;
    QAction *qaAboutAct;

    ProRataMainWidget *prMainWidget;
	QWidget *qwCentralWidget;
	QVBoxLayout *qvblCentralLayout;

    QToolBar *qtbFrequentlyUsedActToolBar;

    QString qsFilename;

	QProgressBar* qpbProgress;
	
	//ProRataSearchDock *prsdSearch;
	QSplashScreen *qssSplash;

};

#endif  //PRORATAMAINWINDOW_H
