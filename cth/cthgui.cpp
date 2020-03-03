/* Main window startup and teardown.
   qtcreator modifies this to add signal/slot handlers.
   These all pass through to the real code in cth_impl.cpp,
   which qtcreator leaves alone.

*/

#include "cthgui.h"
#include "ui_cthgui.h"
#include <iostream>

using namespace std;

CTHgui::CTHgui(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::CTHgui)
{
   ui->setupUi(this);
   setWindowIcon(QIcon(":/cthicon.png"));
   QString title("CTH Clustering App Version: ");
   title = title.append(VERSION);
   setWindowTitle(title);
   setupCtrls();

   QSettings settings("cthgui","cthgui");
   if (settings.contains("geometry"))
   {
      restoreGeometry(settings.value("geometry").toByteArray());
      restoreState(settings.value("windowState").toByteArray());
   }
   timerId = startTimer(1000);  // trigger connect tries
}

CTHgui::~CTHgui()
{
   delete ui;
}



void CTHgui::setupCtrls()
{
   QFont font("courier");

   // button group ids
   ui->linkageButtonGroup->setId(ui->linkageSingle,SINGLE);
   ui->linkageButtonGroup->setId(ui->linkageComplete,COMPLETE);
   ui->linkageButtonGroup->setId(ui->linkageAverage,AVERAGE);
   ui->linkageButtonGroup->setId(ui->linkageWeighted,WEIGHTED);
   ui->linkageButtonGroup->setId(ui->linkageCentroid,CENTROID);
   ui->linkageButtonGroup->setId(ui->linkageMedian,MEDIAN);
   ui->linkageButtonGroup->setId(ui->linkageWard,WARD);
   ui->linkageButtonGroup->setId(ui->linkageArchetype,ARCHETYPE);

   ui->pdistButtonGroup->setId(ui->pdistFromFile,CUSTOM);
   ui->pdistButtonGroup->setId(ui->pdistEuclidean,EUCLIDEAN);
   ui->pdistButtonGroup->setId(ui->pdistSeuclidean,SEUCLIDEAN);
   ui->pdistButtonGroup->setId(ui->pdistCityblock,CITYBLOCK);
   ui->pdistButtonGroup->setId(ui->pdistCosine,COSINE);
   ui->pdistButtonGroup->setId(ui->pdistCorrelation,CORRELATION);
   ui->pdistButtonGroup->setId(ui->pdistSpearman,SPEARMAN);
   ui->pdistButtonGroup->setId(ui->pdistChebychev,CHEBYCHEV);

   ui->basisButtonGroup->setId(ui->basisAuto,1);
   ui->basisButtonGroup->setId(ui->basisPick,2);

   ui->choiceButtonGroup->setId(ui->choiceClick,CLICK);
   ui->choiceButtonGroup->setId(ui->choiceNum,NUM);
   ui->choiceButtonGroup->setId(ui->choiceDist,HDIST);

    // file system browser & CTH file picker
   fsysModel = new QFileSystemModel;
   fsysModel->setRootPath(QDir::currentPath());
   fsysModel->setNameFilterDisables(false);
   fsysModel->setReadOnly(true);
   QStringList filt({"*.cth","*.type"});
   fsysModel->setNameFilters(filt);
   fsysModel->setFilter(QDir::AllDirs | QDir::NoDot | QDir::Files);
   ui->fsysTree->setModel(fsysModel);
   ui->fsysTree->setRootIndex(fsysModel->index(QDir::currentPath()));
   connect(fsysModel, SIGNAL(directoryLoaded(QString)), this,
                             SLOT(model_directoryLoaded(QString)));
   ui->fsysTree->hideColumn(1);
   ui->fsysTree->hideColumn(2);
   ui->fsysTree->hideColumn(3);
   ui->fsysTree->setWindowTitle("Hi");

    // Used to have two different trees, now put all type files in model, drag
    // and drop to destination slots. There may be 3, 4, 5, different types of
    // archtype files, so need more extensible way to manage them.
   centroidModel = new QFileSystemModel;
   centroidModel->setRootPath(QDir::currentPath());
   centroidModel->setNameFilterDisables(false);
   centroidModel->setReadOnly(true);
   QStringList centfilt({"*.type"});
   centroidModel->setNameFilters(centfilt);
   centroidModel->setFilter(QDir::AllDirs | QDir::NoDot | QDir::Files);
   ui->centroidTree->setModel(centroidModel);
   ui->centroidTree->setRootIndex(centroidModel->index(QDir::currentPath()));
   connect(centroidModel, SIGNAL(directoryLoaded(QString)), this,
                            SLOT(centroid_directoryLoaded(QString)));
   ui->centroidTree->hideColumn(1);
   ui->centroidTree->hideColumn(2);
   ui->centroidTree->hideColumn(3);
   ui->centroidTree->setWindowTitle("");
     // these controls are R/O, but we still want the Delete button to work
   QAction* clearAction = ui->archetypeName->findChild<QAction*>();
   if (clearAction)
      clearAction->setEnabled(true); 
   clearAction = ui->swallowName->findChild<QAction*>();
   if (clearAction)
      clearAction->setEnabled(true); 
   clearAction = ui->lareflexName->findChild<QAction*>();
  if (clearAction)
      clearAction->setEnabled(true); 

   if (!debug)
      ui->guiConnect->hide();  // hide this from most users
   ui->guiInfo->setReadOnly(true);
   ui->guiInfo->setFont(font);
}

// File menu quit handler
void CTHgui::on_actionClose_triggered()
{
   doQuit();
}

void CTHgui::on_basisAuto_clicked()
{
   basisAuto();
}

void CTHgui::on_basisPick_clicked()
{
    basisPick();
}

void CTHgui::on_guiConnect_clicked()
{
    guiConnect();
}

void CTHgui::timerEvent(QTimerEvent *event)
{
   tryConnect(event->timerId());
}

void CTHgui::on_sendOctave_clicked()
{
   sendControls();
}

void CTHgui::on_choiceClick_clicked()
{
   choiceClick();
}

void CTHgui::on_choiceNum_clicked()
{
   choiceNum();
}

void CTHgui::on_sendCmd_clicked()
{
   sendCmd();
}

void CTHgui::on_doQuit_clicked()
{
    doQuit();
}

void CTHgui::closeEvent(QCloseEvent *evt)
{
    doCloseQuit(evt);
}

// str is absolute path to .. 
void CTHgui::model_directoryLoaded(QString str)
{
   dirLoaded(str);
}

void CTHgui::centroid_directoryLoaded(QString str)
{
   centroidDirLoaded(str);
}

void CTHgui::on_filePath_returnPressed()
{
    filePathReturn();
}

void CTHgui::on_fsysTree_activated(const QModelIndex &index)
{
   fsysTreeSel(index);
}

void CTHgui::on_fsysTree_doubleClicked(const QModelIndex &index)
{
   fsysTreeSel(index);
}

void CTHgui::on_choiceDist_clicked()
{
   choiceDist();
}

void CTHgui::on_sendParam_clicked()
{
   sendParam();
}

void CTHgui::on_inconsistencyCalc_toggled(bool checked)
{
   inconsistencyCalc(checked);
}

void CTHgui::on_selectStereo_clicked()
{
   selectStereo();
}

void CTHgui::on_saveClusters_clicked()
{
   saveClusters();
}

void CTHgui::on_startOver_clicked()
{
   startOver();
}

void CTHgui::on_sendCmdText_returnPressed()
{
   on_sendCmd_clicked();
}

void CTHgui::on_scrollFwd_clicked()
{
   scrollFwd();
}

void CTHgui::on_scrollBack_clicked()
{
   scrollBack();
}

void CTHgui::on_selAll_clicked()
{
   selAll();
}

void CTHgui::on_selNone_clicked()
{
   selNone();
}

void CTHgui::on_getCTHInfo_clicked()
{
   getCTHInfo();
}

void CTHgui::on_basisPCA_clicked()
{
    basisPCA();
}

void CTHgui::on_plot2D_clicked(bool checked)
{
    plot2D(checked);
}

void CTHgui::on_exportAtlas_clicked()
{
    exportAtlas();
}

void CTHgui::on_getCTHInfoText_returnPressed()
{
   on_getCTHInfo_clicked();
}

void CTHgui::on_saveClustersData_returnPressed()
{
   on_saveClusters_clicked();
}

void CTHgui::on_sendparamText_returnPressed()
{
    on_sendParam_clicked();
}

void CTHgui::on_selectStereoData_returnPressed()
{
   on_selectStereo_clicked();
}


void CTHgui::on_exportAtlasText_returnPressed()
{
    exportAtlas();
}

void CTHgui::on_fuzzyCluster_clicked(bool checked)
{
   fuzzyClick(checked);
}

void CTHgui::on_exportCentroids_clicked()
{
   exportCentroids();
}

void CTHgui::on_exportCentroidsText_returnPressed()
{
   exportCentroids();
}

void CTHgui::on_actionHelp_triggered()
{
   doHelp();
}

void CTHgui::on_actionQuit_triggered()
{
   doQuit();
}

void CTHgui::on_centroidTree_activated(const QModelIndex &index)
{
   centroidTreeSel(index);
}

void CTHgui::on_centroidTree_doubleClicked(const QModelIndex &index)
{
   centroidTreeSel(index);
}

void CTHgui::on_centroidPath_returnPressed()
{
   centroidFilePathReturn();
}


void CTHgui::on_closeInfo_clicked()
{
   closeInfo();
}

void CTHgui::on_exportPairs_clicked()
{
   exportPairs();
}
