/* Implementation of UI controls

*/

#include <QObject>
#include <QMessageBox>
#include <QRadioButton>
#include <QListWidget>
#include <qsocketnotifier.h>

#include <QGuiApplication>
#include <QWindow>
#include <QScreen>

#include <stdio.h>
#include <iostream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>

#include "cthgui.h"
#include "helpbox.h"
#include "ui_cthgui.h"

#pragma GCC diagnostic ignored "-Wunused-result"

using namespace std;

const QString END("_END_");
const QString input_req("INPUT REQUIRED");
const QString param_req("PARAM REQUIRED");
const QString error_warn("ERROR WARN");

// Send our current selections to the octave server waiting to hear from us
void CTHgui::sendOctave(SENDTYPE send)
{
   QString cmd;

   switch (send)
   {
      case CMD:  // send a command
         cmdStr(cmd);
         break;
      case PARAMS:
         paramStr(cmd);
         break;
      case STEREO_CLUSTERS:
         stereoStr(cmd);
         break;
      case SAVE_CLUSTERS:
         saveStr(cmd);
         break;
      case START_OVER:
         restartStr(cmd);
         break;
      case FORWARD:
         forwardStr(cmd);
         break;
      case BACKWARD:
         backwardStr(cmd);
         break;
      case CTHINFO:
         cthinfoStr(cmd);
         break;
      case EXPORT_ATLAS:
         exportStr(cmd);
         break;
      case EXPORT_ARCHETYPES:
         exportCents(cmd);
         break;
      case CONTROLS:
         controlsStr(cmd);
         break;
      case CLOSE_INFO:
         closeInfoFigs(cmd);
         break;
      case EXPORT_PAIRS:
         exportPairsToFile(cmd);
         break;
      default:
         break;
   }

   if (!cmd.isEmpty())
   {
      ui->guiInfo->repaint();
      QByteArray to_send = cmd.toLatin1();  // unicode to latin, this probably
                                            // fails in some locales
      int result = ::send(octaveSocket,to_send,to_send.length(),0);
      if (result != (int) cmd.length())
        uiWarn("WARNING:  All data not sent to octave.\n");
   }
}


void CTHgui::cmdStr(QString& cmd)
{
   QString curr = ui->sendCmdText->text();
   QTextStream(&cmd) << "CMD " << curr << " END";
}

void CTHgui::paramStr(QString& cmd)
{
   QString param = ui->sendparamText->text();
   if (param.length() == 0)
      uiWarn("No parameter to send, try again.\n");
   else
   {
      QTextStream(&cmd) << "PARAMS " << param << " END";
      ui->sendOctave->setEnabled(true);  // can use these again
      ui->sendCmd->setEnabled(true);
      ui->selectStereo->setEnabled(true);
      ui->saveClusters->setEnabled(true);
      ui->startOver->setEnabled(true);
      ui->scrollFwd->setEnabled(true);
      ui->scrollBack->setEnabled(true);
      ui->sendParam->setEnabled(false);
      ui->getCTHInfo->setEnabled(true);
      ui->exportAtlas->setEnabled(true);
      ui->exportPairs->setEnabled(true);
      ui->closeInfo->setEnabled(haveCthInfo);
   }
}

void CTHgui::stereoStr(QString& cmd)
{
   QString clusts = ui->selectStereoData->text();
   if (clusts.isEmpty())
   {
      uiWarn("No clusters to send, try again.\n");
   }
   else
      QTextStream(&cmd) << "STEREO_CLUSTERS " << clusts << " END";
}

void CTHgui::saveStr(QString& cmd)
{
   QString clusts = ui->saveClustersData->text();
   if (clusts.isEmpty())
      uiWarn("No clusters to send, try again.\n");
   else
      QTextStream(&cmd) << "SAVE_CLUSTERS " << clusts << " END";
}

void CTHgui::restartStr(QString& cmd)
{
   QTextStream(&cmd) << "START_OVER END";
   ui->selectStereo->setEnabled(false);
   ui->saveClusters->setEnabled(false);
   ui->startOver->setEnabled(false);
   ui->scrollFwd->setEnabled(false);
   ui->scrollBack->setEnabled(false);
   ui->getCTHInfo->setEnabled(false);
   ui->exportAtlas->setEnabled(false);
   ui->exportCentroids->setEnabled(false);
   ui->exportPairs->setEnabled(false);
   haveCthInfo = false;
   ui->closeInfo->setEnabled(haveCthInfo);
   ui->getCTHInfoText->clear();
   ui->saveClustersData->clear();
   ui->sendparamText->clear();
   ui->selectStereoData->clear();
   ui->exportAtlasText->clear();
   ui->exportPairsText->clear();
   ui->exportCentroidsText->clear();
   ui->sendCmdText->clear();
}

void CTHgui::forwardStr(QString& cmd)
{
   QTextStream(&cmd) << "FORWARD END";
}

void CTHgui::backwardStr(QString& cmd)
{
   QTextStream(&cmd) << "BACKWARD END";
}

void CTHgui::cthinfoStr(QString& cmd)
{
   QString info = ui->getCTHInfoText->text();
   if (info.isEmpty())
      uiWarn("No CTH numbers to send, try again.\n");
   else
   {
      QTextStream(&cmd) << "CTHINFO " << info << " END";
      haveCthInfo = false;
      ui->closeInfo->setEnabled(haveCthInfo);
   }
}

void CTHgui::exportStr(QString& cmd)
{
   QString fname = ui->exportAtlasText->text();
   if (fname.isEmpty())
      uiWarn("No filename to save to, try again.\n");
   else
      QTextStream(&cmd) << "EXPORT_ATLAS " << fname << " END";
}

void CTHgui::exportCents(QString& cmd)
{
   QString fname = ui->exportCentroidsText->text();
   if (fname.isEmpty())
      uiWarn("No filename to save to, try again.\n");
   else
      QTextStream(&cmd) << "EXPORT_ARCHETYPES " << fname << " END";
}

void CTHgui::exportPairsToFile(QString& cmd)
{
   QString fname = ui->exportPairsText->text();
   if (fname.isEmpty())
      uiWarn("No filename to save to, try again.\n");
   else
      QTextStream(&cmd) << "EXPORT_PAIRS " << fname << " END";
}


void CTHgui::closeInfoFigs(QString& cmd)
{
   QTextStream(&cmd) << "CLOSE_INFO END";
}

// This gets called when the Create Plots or the Cluster Using Archetypes 
// button is clicked.
void CTHgui::controlsStr(QString& cmd)
{
   bool have_errs = false;
   bool ok;
   QString num;
   CLUSTS pickclust;
   int numclust = 0;
   BASIS basis;
   int p0, p1, p2, p3;
   char *fname = nullptr;
   char *type_name = nullptr;
   char *swall_name = nullptr;
   char *lareflex_name = nullptr;
   float dens_percent = 2.0;
   float fuzzy_cut = 0.1;
   int inconsistent;
   int incon_depth=2;
   int link_algo;
   p0 = p1 = p2 = p3 = 0;
   pid_t pid = getpid();
   QString file;
   WId winid=winId();

    // do items that can be error-checked first
   num = ui->fuzzyClusterCutoff->text();
   fuzzy_cut = num.toFloat(&ok);
   if (!ok || fuzzy_cut <= 0.0 || fuzzy_cut >= 1.0) 
      fuzzy_cut = 0;

   inconsistent = ui->inconsistencyCalc->isChecked();
   if (inconsistent)
   {
      num  = ui->inconsistentDepth->text();
      incon_depth = num.toInt(&ok);
      if (!ok || incon_depth < 0)
      {
         uiWarn("Inconsistency cutoff is not valid.");
         have_errs = true;
      }
   }

   pickclust = CLUSTS(ui->choiceButtonGroup->checkedId());
   switch (pickclust)
   {
      case NUM:
         num = ui->choiceEnterNum->text();
         numclust = num.toInt(&ok);
           // any non-valid input is assumed to indicate that
           // octave should prompt for a value later.
         if (!ok || numclust < 0)
         {
            ui->choiceEnterNum->setText("0");
            numclust=0;
         }
         break;
      case HDIST:
         num = ui->choiceEnterDist->text();
         numclust = num.toInt(&ok);
         if (!ok || numclust < 0)
         {
            ui->choiceEnterDist->setText("0");
            numclust=0;
         }
         break;

      case CLICK:
      default:
         break;
   }

   if (ui->basisAuto->isChecked())
      basis = BASIS_MOSTDIST;
   else if (ui->basisPCA->isChecked())
      basis = BASIS_PCA;
   else
   {
      basis = BASIS_PICK;
      num  = ui->coordP0->text();
      p0 = num.toInt(&ok);
      if (!ok || p0 < 1)
      {
         uiWarn("Point p0 is not valid.");
         have_errs = true;
      }
      num  = ui->coordP1->text();
      p1 = num.toInt(&ok);
      if (!ok || p1 < 1)
      {
         uiWarn("Point p1 is not valid.");
         have_errs = true;
      }
      num  = ui->coordP2->text();
      p2 = num.toInt(&ok);
      if (!ok || p2 < 1)
      {
         uiWarn("Point p2 is not valid.");
         have_errs = true;
      }
      if (!ui->plot2D->isChecked())
      { 
         num  = ui->coordP3->text();
         p3 = num.toInt(&ok);
         if (!ok || p3 < 1)
         {
            uiWarn("Point p3 is not valid.");
            have_errs = true;
         }
      }
   }

   QByteArray array;
   QModelIndex currindex = ui->fsysTree->currentIndex();
   QFileInfo info;
   if (!currindex.isValid())
   {
      uiWarn("Select a file to load");
      have_errs = true; 
   }
   else
   {
      info = fsysModel->fileInfo(currindex);
      bool dir = info.isDir();
      if (dir)
      {
         uiWarn("The selected item is a directory");
         have_errs = true; 
      }
      else
      {
         file = info.absoluteFilePath();
         file.replace(" ","<<**SPACE**>>");   // meta-info because spaces break everything
         array = file.toLocal8Bit();
         asprintf(&fname,"%s", array.data());
            // create default export csv fname
         ui->exportAtlasText->setText(info.completeBaseName());
         ui->exportCentroidsText->setText(info.completeBaseName());
         ui->exportPairsText->setText(info.completeBaseName());
      }
   }

   link_algo = ui->linkageButtonGroup->checkedId();
   if (link_algo == ARCHETYPE)
   {
      QByteArray array;
      QFileInfo info;
      QString file;
      QString path = ui->centroidPath->text();

       // If using archetypes for clustering, we don't know what .type files we
       // need because (so far) we never open the file to look at it, so we
       // can't check for errors here. We don't even know if the .type files
       // are appropriate. Forward all this to the octave code and let it sort
       // it out.
         // standard archetypes
      file = ui->archetypeName->text();
      if (file.length() == 0)
         asprintf(&type_name,"%s", "NONE");
      else
      {
         file = path + "/" + ui->archetypeName->text();
         file.replace(" ","<<**SPACE**>>");   // meta-info because spaces break everything
         array = file.toLocal8Bit();
         asprintf(&type_name,"%s", array.data());
      }
         // swallow archetypes
      asprintf(&swall_name,"%s", "NONE");  // assume the worst
      file = ui->swallowName->text();
      if (file.length() > 0)
      {
         file = path + "/" + ui->swallowName->text();
         file.replace(" ","<<**SPACE**>>");
         array = file.toLocal8Bit();
         asprintf(&swall_name,"%s", array.data());
      }
         // lareflex archetypes
      asprintf(&lareflex_name,"%s", "NONE");  // assume the worst
      file = ui->lareflexName->text();
      if (file.length() > 0)
      {
         file = path + "/" + ui->lareflexName->text();
         file.replace(" ","<<**SPACE**>>");
         array = file.toLocal8Bit();
         asprintf(&lareflex_name,"%s", array.data());
      }
   }
   else
   {
      asprintf(&type_name,"NONE");
      asprintf(&swall_name,"NONE");
      asprintf(&lareflex_name,"NONE");
   }
   if (link_algo == ARCHETYPE && ui->pdistButtonGroup->checkedId() == CUSTOM)
   {
      ui->pdistEuclidean->setChecked(true);
      uiWarn("The Custom distance metric is not available for archetypes.\nDefaulting to Euclidean distance.");
   }
   int pd_algo = ui->pdistButtonGroup->checkedId();
   int plotdim = ui->plot2D->isChecked() ? DIM2 : DIM3;
   int cophen = ui->copheneticCalc->isChecked();
   int errorbars = ui->errorbarsShow->isChecked();
   int curveshow = ui->curveShow->isChecked();
   int curvederive = ui->curveDerive->isChecked();
   int errscale = ui->errPlotScale->isChecked();
   int kmeans = ui->kmeansPlot->isChecked();
   int autoscale = ui->varBarY->isChecked();
   int flatplots = ui->flatPlots->isChecked();
   int stereoper = ui->stereoPerExp->isChecked();
   int pairwise = ui->pairWise->isChecked();
   int typepairs = ui->typePairs->isChecked();
   int log_dist = ui->logDist->isChecked();
   int isoplots = ui->isoMap->isChecked();
   int isoover = ui->isoMapOverlay->isChecked();
   int isoembed = ui->isoMapEmbed->isChecked();
   int isocolor = ui->isoMapColor->isChecked();
   int isoscale = ui->isoVarBarY->isChecked();
   int mod1 = ui->modDepth1->isChecked();
   int mod2 = ui->modDepth2->isChecked();
   int fuzzy = ui->fuzzyCluster->isChecked();
   int scrn_cap = ui->screenCap->isChecked();
   int no_flats = ui->noFlats->isChecked();
   int color;
   if (ui->colorOn->isChecked())
      color = 2;
   else
      color = 1;

   if (!have_errs)
   {
      QTextStream(&cmd) << "CONTROLS " << "0"
      << " fname " << fname 
      << " typename " << type_name 
      << " pd_algo " << pd_algo 
      << " link_algo " << link_algo 
      << " color " << color 
      << " densperc " << dens_percent 
      << " cophen " << cophen 
      << " incon_sel " << inconsistent 
      << " incon_depth " << incon_depth 
      << " errorbars " << errorbars 
      << " errscale " << errscale 
      << " kmeans " << kmeans 
      << " plotdim " << plotdim 
      << " basis " << basis 
      << " p0 " << p0 
      << " p1 " << p1 
      << " p2 " << p2 
      << " p3 " << p3 
      << " pickclust " << pickclust 
      << " numclust " << numclust 
      << " pid " << pid 
      << " flatplots " << flatplots 
      << " stereoper " << stereoper 
      << " pairwise " << pairwise 
      << " typepairs " << typepairs 
      << " log_dist " << log_dist 
      << " autoscale " << autoscale 
      << " winid " << winid 
      << " isoplots " << isoplots 
      << " isoover " << isoover 
      << " isoembed " << isoembed 
      << " isocolor " << isocolor 
      << " isoscale " << isoscale 
      << " mod1 " << mod1 
      << " mod2 " << mod2 
      << " curveshow " << curveshow 
      << " curvederive " << curvederive 
      << " fuzzy " << fuzzy
      << " fuzzycutoff " << fuzzy_cut
      << " scrncap " << scrn_cap
      << " swallname " << swall_name 
      << " lareflexname " << lareflex_name 
      << " noflats " << no_flats 
      << " end END"; 
      ui->selectStereo->setEnabled(true); // need some plots for these to make sense
      ui->saveClusters->setEnabled(true);
      ui->startOver->setEnabled(true);
      ui->scrollFwd->setEnabled(true);
      ui->scrollBack->setEnabled(true);
      ui->getCTHInfo->setEnabled(true);
      ui->exportAtlas->setEnabled(true);
      ui->exportCentroids->setEnabled(true);
      ui->exportPairs->setEnabled(true);
      ui->closeInfo->setEnabled(haveCthInfo);
   }
   free(fname);
   free(type_name);
}


void CTHgui::sendControls()
{
   sendOctave(CONTROLS);
}

void CTHgui::sendCmd()
{
   sendOctave(CMD);
}

void CTHgui::sendParam()
{
   sendOctave(PARAMS);
}

void CTHgui::selectStereo()
{
   sendOctave(STEREO_CLUSTERS);
}

void CTHgui::saveClusters()
{
   sendOctave(SAVE_CLUSTERS);
}

void CTHgui::startOver()
{
   sendOctave(START_OVER);
}

void CTHgui::scrollFwd()
{
   sendOctave(FORWARD);
}

void CTHgui::scrollBack()
{
   sendOctave(BACKWARD);
}

void CTHgui::getCTHInfo()
{
   sendOctave(CTHINFO);
   haveCthInfo = true;
   ui->closeInfo->setEnabled(haveCthInfo);
}

void CTHgui::exportAtlas()
{
   sendOctave(EXPORT_ATLAS);
}

void CTHgui::closeInfo()
{
   sendOctave(CLOSE_INFO);
   haveCthInfo = false;
   ui->closeInfo->setEnabled(haveCthInfo);
}

void CTHgui::exportCentroids()
{
   sendOctave(EXPORT_ARCHETYPES);
}

void CTHgui::exportPairs()
{
   sendOctave(EXPORT_PAIRS);
}

void CTHgui::basisAuto()
{
   ui->coordP0->setEnabled(false);
   ui->coordP1->setEnabled(false);
   ui->coordP2->setEnabled(false);
   ui->coordP3->setEnabled(false);
}

void CTHgui::basisPCA()
{
   ui->coordP0->setEnabled(false);
   ui->coordP1->setEnabled(false);
   ui->coordP2->setEnabled(false);
   ui->coordP3->setEnabled(false);
}

void CTHgui::basisPick()
{
   ui->coordP0->setEnabled(true);
   ui->coordP1->setEnabled(true);
   ui->coordP2->setEnabled(true);
   if (ui->plot2D->isChecked() == true)
      ui->coordP3->setEnabled(false);
   else
      ui->coordP3->setEnabled(true);
}


void CTHgui::inconsistencyCalc(bool checked)
{
   if (checked)
      ui->inconsistentDepth->setEnabled(true);
   else
      ui->inconsistentDepth->setEnabled(false);
}

void CTHgui::plot2D(bool checked)
{
   if (checked)
   {
      ui->coordP3->setEnabled(false);
      ui->coordP3->setText("");
      ui->basisPick->setText("Pick 3 Points");
   }
   else
   {
      if (ui->basisPick->isChecked())
         ui->coordP3->setEnabled(true);
      ui->basisPick->setText("Pick 4 Points");
   }
}

void CTHgui::choiceClick()
{
   ui->choiceEnterNum->setEnabled(false);
   ui->choiceEnterDist->setEnabled(false);
}

void CTHgui::choiceNum()
{
   ui->choiceEnterNum->setEnabled(true);
   ui->choiceEnterNum->setFocus();
   ui->choiceEnterDist->setEnabled(false);
}

void CTHgui::choiceDist()
{
   ui->choiceEnterDist->setEnabled(true);
   ui->choiceEnterNum->setEnabled(false);
}

void CTHgui::fuzzyClick(bool checked)
{
   ui->fuzzyClusterCutoff->setEnabled(checked);
}

void CTHgui::selAll()
{
   ui->errorbarsShow->setChecked(true);
   ui->flatPlots->setChecked(true);
   ui->isoMap->setChecked(true);
   ui->isoMapEmbed->setChecked(true);
   ui->isoMapColor->setChecked(true);
   ui->isoMapOverlay->setChecked(true);
   ui->kmeansPlot->setChecked(true);
   ui->logDist->setChecked(true);
   ui->pairWise->setChecked(true);
   ui->typePairs->setChecked(true);
   ui->curveShow->setChecked(true);
   ui->curveDerive->setChecked(true);
   ui->fuzzyCluster->setChecked(true);
   ui->fuzzyClusterCutoff->setEnabled(true);
}

void CTHgui::selNone()
{
   ui->errorbarsShow->setChecked(false);
   ui->flatPlots->setChecked(false);
   ui->isoMap->setChecked(false);
   ui->isoMapEmbed->setChecked(false);
   ui->isoMapColor->setChecked(false);
   ui->isoMapOverlay->setChecked(false);
   ui->kmeansPlot->setChecked(false);
   ui->logDist->setChecked(false);
   ui->pairWise->setChecked(false);
   ui->typePairs->setChecked(false);
   ui->curveShow->setChecked(false);
   ui->curveDerive->setChecked(false);
   ui->fuzzyCluster->setChecked(false);
   ui->fuzzyClusterCutoff->setEnabled(false);
}

void CTHgui::guiConnect()
{
   tryConnect(0);
}

// This can be called when a timer event occurs, or as
// a result of a button click.
void CTHgui::tryConnect(int timerId)
{
   int sock_fd;
   struct sockaddr_in serv_addr;
   struct hostent *server;
   static int give_up = 5;
   char msg[1024];

   sock_fd = socket(AF_INET,SOCK_STREAM,0);
   if (sock_fd < 0)
   {
      uiWarn("error opening socket");
      return;
   }

   server = gethostbyname("localhost");
   if (server == nullptr)
   {
      uiWarn("Can't find localhost");
      return;
   }

   bzero((char*) &serv_addr,sizeof(serv_addr));
   serv_addr.sin_family = AF_INET;
   bcopy((char *)server->h_addr,(char*) &serv_addr.sin_addr.s_addr,server->h_length);
   serv_addr.sin_port = htons(octavePort);
   sprintf(msg,"Trying to connect to port %d. . .",octavePort);
   ui->guiInfo->append(msg);
   repaint();
   ui->guiInfo->repaint();
   if (::connect(sock_fd,(struct sockaddr*)&serv_addr,sizeof(serv_addr)) < 0)
   {
      uiWarn("Connect failed. . .");
      ui->guiInfo->repaint();
      ::close(sock_fd);
      --give_up;
      if (!give_up)
      {
         uiWarn("Connect failed, giving up.");
         ui->guiInfo->append("If debugging, use button to connect manually.");
         ui->guiConnect->show();
         ui->guiConnect->setEnabled(true);
         if (timerId)
         {
            killTimer(timerId);
            timerId = 0;
         }
         give_up = 5;
         ui->guiInfo->repaint();
      }
   }
   else
   {
      ui->guiInfo->append("Connected!\n");
      octaveSocket = sock_fd;
      if (timerId)
      {
         killTimer(timerId);
         timerId = 0;
      }
      ui->sendOctave->setEnabled(true);
      ui->sendCmd->setEnabled(true);
      ui->doQuit->setEnabled(true);
      ui->guiConnect->setEnabled(false);
      ui->guiConnect->hide();

      QSocketNotifier *sn;    // set up notifier callback
      sn = new QSocketNotifier(sock_fd,QSocketNotifier::Read);
      QObject::connect(sn,SIGNAL(activated(int)),this,SLOT(dataReceived()) );
   }
}

/* This gets called when the octave app sends us stuff to print, or other
   commands.  
   Sometimes the msgs get fragmented, so we expect each msg to be terminate
   with a special marker (at this time "__END__").  Obviously if this is part
   of a message, things break.  We accumulate text until we see the end marker.
   The append() method will add a newline.
*/

static QString curr_oct_msg;

void CTHgui::dataReceived()
{
   char oct_msg_frag[4096];
   int prefix;
   char *str = nullptr;
   int str_end;
   QString oct_msg;

   int count = ::recv(octaveSocket,oct_msg_frag,sizeof(oct_msg_frag),MSG_DONTWAIT);
   if (count > 0)
   {
      oct_msg_frag[count]=0;  // just bytes, terminate so it's a string
      curr_oct_msg.append(oct_msg_frag);
      str_end = curr_oct_msg.indexOf(END);
      if (str_end == -1)  // no end yet, continue to accumulate
         return;
      else
      {
         oct_msg = curr_oct_msg.mid(0,str_end-1); 
         curr_oct_msg.clear();      
      }
   }
   else
      return;


    // we differentiate between info and requests for input 
    // with a fixed prefix string, but we don't print it.
   if ((prefix = oct_msg.indexOf(input_req)) != -1)
   {
      QColor color = ui->guiInfo->textColor();
      int wght = ui->guiInfo->fontWeight();
      ui->guiInfo->setTextColor(Qt::blue);
      ui->guiInfo->setFontWeight(QFont::Bold);
      ui->guiInfo->append(oct_msg.mid(input_req.length(),-1));
      ui->guiInfo->setTextColor(color);
      ui->guiInfo->setFontWeight(wght);
   }
   else if ((prefix = oct_msg.indexOf(param_req)) != -1)
   {
      ui->sendParam->setEnabled(true);
      ui->sendOctave->setEnabled(false);
      ui->sendCmd->setEnabled(false);
      ui->selectStereo->setEnabled(false);
      ui->saveClusters->setEnabled(false);
      QColor color = ui->guiInfo->textColor();
      int wght = ui->guiInfo->fontWeight();
      ui->guiInfo->setTextColor(Qt::blue);
      ui->guiInfo->setFontWeight(QFont::Bold);
      QString param_text = oct_msg.mid(param_req.length(),-1);
      QByteArray c_stng = param_text.toLatin1();
      asprintf(&str,"Enter parameter(s) for%s in the Send Param(s) edit box to the left.", c_stng.data()); 
      ui->guiInfo->append(str);
      free(str);
      ui->guiInfo->append("Then click on Send Param.");
      ui->guiInfo->setTextColor(color);
      ui->guiInfo->setFontWeight(wght);
   }
   else if ((prefix = oct_msg.indexOf(error_warn)) != -1)
   {
      uiWarn(oct_msg.mid(error_warn.length(),-1));
   }
   else
      ui->guiInfo->append(oct_msg);
}

void CTHgui::doQuit()
{
   close();
}

/* The Wrap Up function.
   This sends a message to the octave side, saves some
   window settings, and eventually the app exits.
*/
void CTHgui::doCloseQuit(QCloseEvent* evt)
{
   char str[1024];
   pid_t pid = getpid();

   sprintf(str,"QUIT %d END",pid);
   ::send(octaveSocket,str,strlen(str),0);  // (do we need to delay here?)

   QSettings settings("cthgui","cthgui");
   settings.setValue("geometry",saveGeometry());
   settings.setValue("windowState",saveState());
   evt->accept();
}

// User has made selection in fsys browser, update the contents
void CTHgui::fsysTreeSel(const QModelIndex &index)
{
   QString filePath = fsysModel->fileInfo(index).absoluteFilePath();
   ui->fsysTree->setRootIndex(fsysModel->setRootPath(filePath));
}

// user types in a path. ENTER triggers this. If valid, switch to it
void CTHgui::filePathReturn()
{
   QString fpath = ui->filePath->text();
   if (QDir(fpath).exists())
      ui->fsysTree->setRootIndex(fsysModel->setRootPath(fpath));
}

// dir has changed, update path box
void CTHgui::dirLoaded(QString str)
{
   ui->filePath->setText(str);
}

// User has made selection in archetype browser, update the contents
void CTHgui::centroidTreeSel(const QModelIndex &index)
{
   QString filePath = centroidModel->fileInfo(index).absoluteFilePath();
   ui->centroidTree->setRootIndex(centroidModel->setRootPath(filePath));
}

void CTHgui::centroidFilePathReturn()
{
   QString fpath = ui->centroidPath->text();
   if (QDir(fpath).exists())
      ui->centroidTree->setRootIndex(centroidModel->setRootPath(fpath));
}

// Pick a default Archetype .type file if we see one we like
// and update current path control.
void CTHgui::centroidDirLoaded(QString /* str */)
{
   QApplication::sendPostedEvents(); // control may not be completely inited, wait for it

   QString path = centroidModel->rootPath();
   ui->centroidPath->setText(path);
   QModelIndex root = centroidModel->index(path);
   int rows=centroidModel->rowCount(root);
   QModelIndex idx;
   for (int row=0;row<rows;++row)
   {
      idx= centroidModel->index(row,0,root);
      if(idx.isValid())
      {
         QString name = idx.data(QFileSystemModel::FileNameRole).toString();
         if (name.contains(QRegularExpression("^Arch.*.type$")))
         {
            ui->centroidTree->setCurrentIndex(idx);
            ui->centroidTree->scrollTo(idx,QAbstractItemView::PositionAtTop);
            ui->centroidTree->expand(idx);
            break;
         }
      }
   }
}

void CTHgui::doHelp()
{
   helpbox *help = new helpbox(this);
   help->exec();
}

void CTHgui::uiWarn(QString msg)
{
   QColor color = ui->guiInfo->textColor();
   int wght = ui->guiInfo->fontWeight();
   ui->guiInfo->setTextColor(Qt::red);
   ui->guiInfo->setFontWeight(QFont::Bold);
   ui->guiInfo->append(msg);
   ui->guiInfo->setTextColor(color);
   ui->guiInfo->setFontWeight(wght);
}


