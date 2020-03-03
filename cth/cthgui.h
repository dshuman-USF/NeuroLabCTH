#ifndef CTHGUI_H
#define CTHGUI_H

#include <QMainWindow>
#include <QDir>
#include <QModelIndex>
#include <QFileSystemModel>
#include <QCloseEvent>
#include <QTextStream>
#include <QSettings>
#include <getopt.h>

namespace Ui {
class CTHgui;
}

class CTHgui : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit CTHgui(QWidget *parent = 0);
    virtual ~CTHgui();

   // these are ordered in the octave function parm order, or in
   // the order cth_project.m or other scripts expect.
   enum LINKAGE {SINGLE=1, COMPLETE, AVERAGE, WEIGHTED, CENTROID, MEDIAN, WARD, ARCHETYPE};
     // note:  There are some dist functions that are not useful for us.
     // The mahalanobis function in octave seems very broken. 
     // Also, HAMMING and JACCARD do nothing for us, they seem
     // to be meant for very different kinds of data, such as text.  Minkowski
     // is redundant since we don't allow input of P, so it is always euclidean.
   enum PDIST {CUSTOM=1, EUCLIDEAN, SEUCLIDEAN, CITYBLOCK, 
               COSINE, CORRELATION, SPEARMAN, CHEBYCHEV};
   enum PLOTDIM {DIM2=1, DIM3};
   enum CLUSTS {NUM=1, HDIST, CLICK};
   enum SENDTYPE {CONTROLS, CMD, PARAMS, STEREO_CLUSTERS,SAVE_CLUSTERS,START_OVER,FORWARD,BACKWARD,CTHINFO,EXPORT_ATLAS,EXPORT_ARCHETYPES,CLOSE_INFO,EXPORT_PAIRS};
   enum BASIS {BASIS_MOSTDIST=1,BASIS_PCA,BASIS_PICK};

   void setPort(int port) { octavePort = port;};
   void setDebug(bool db) { debug = db; };

private slots:
    void on_actionClose_triggered();
    void on_basisAuto_clicked();
    void on_basisPick_clicked();
    void on_plot2D_clicked(bool);
    void on_guiConnect_clicked();
    void on_sendOctave_clicked();
    void dataReceived();
    void on_choiceClick_clicked();
    void on_choiceNum_clicked();
    void on_sendCmd_clicked();
    void on_doQuit_clicked();
    void model_directoryLoaded(QString);
    void centroid_directoryLoaded(QString);
    void on_filePath_returnPressed();
    void on_fsysTree_activated(const QModelIndex &index);
    void on_fsysTree_doubleClicked(const QModelIndex &index);
    void on_choiceDist_clicked();
    void on_sendParam_clicked();
    void on_inconsistencyCalc_toggled(bool checked);
    void on_selectStereo_clicked();
    void on_saveClusters_clicked();
    void on_startOver_clicked();
    void on_sendCmdText_returnPressed();
    void on_scrollFwd_clicked();
    void on_scrollBack_clicked();
    void on_selAll_clicked();
    void on_selNone_clicked();
    void on_getCTHInfo_clicked();
    void on_basisPCA_clicked();
    void on_getCTHInfoText_returnPressed();
    void on_saveClustersData_returnPressed();
    void on_sendparamText_returnPressed();
    void on_selectStereoData_returnPressed();
    void on_exportAtlas_clicked();
    void on_exportAtlasText_returnPressed();
    void on_fuzzyCluster_clicked(bool checked);
    void on_exportCentroids_clicked();
    void on_exportCentroidsText_returnPressed();
    void on_actionHelp_triggered();
    void on_actionQuit_triggered();
    void on_centroidTree_activated(const QModelIndex &index);
    void on_centroidTree_doubleClicked(const QModelIndex &index);
    void on_centroidPath_returnPressed();
    void on_closeInfo_clicked();
    void on_exportPairs_clicked();

protected:
    void timerEvent(QTimerEvent *event);
    void closeEvent(QCloseEvent *evt);

private:
    Ui::CTHgui *ui;

   int timerId;
   int octavePort;
   int octaveSocket;
   bool debug;
   bool haveCthInfo=false;

   QFileSystemModel *fsysModel;
   QFileSystemModel *centroidModel;

   void setupCtrls();

     // impl of signal handlers
   void sendOctave(SENDTYPE);
   void basisAuto();
   void basisPCA();
   void basisPick();
   void plot2D(bool);
   void guiConnect();
   void tryConnect(int);
   void choiceClick();
   void fuzzyClick(bool);
   void choiceNum();
   void choiceDist();
   void doQuit();
   void doCloseQuit(QCloseEvent* );
   void fsysTreeSel(const QModelIndex &index);
   void filePathReturn();
   void dirLoaded(QString);
   void centroidFilePathReturn();
   void centroidDirLoaded(QString);
   void centroidTreeSel(const QModelIndex &index);
   void inconsistencyCalc(bool);
   void sendParam();
   void sendCmd();
   void sendControls();
   void selectStereo();
   void saveClusters();
   void startOver();
   void scrollFwd();
   void scrollBack();
   void selAll();
   void selNone();
   void getCTHInfo();
   void exportAtlas();
   void exportCentroids();
   void exportPairs();
   void closeInfo();
   void cmdStr(QString& );
   void paramStr(QString& );
   void stereoStr(QString& );
   void saveStr(QString& );
   void restartStr(QString& );
   void forwardStr(QString& );
   void backwardStr(QString& );
   void cthinfoStr(QString& );
   void exportStr(QString& );
   void controlsStr(QString&);
   void exportCents(QString& );
   void closeInfoFigs(QString& );
   void exportPairsToFile(QString& );
   void doHelp();
   void uiWarn(QString msg);
};

#endif // CTHGUI_H
