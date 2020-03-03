// This draws a maximized window on each physical monitor/screen and
// obtains frame dimensions, the usable area on the screen, and screen res.

#include <QApplication>
#include <QGuiApplication>
#include <QWindow>
#include <QScreen>
#include <QtWidgets/QSystemTrayIcon>
#include <QtWidgets/QDesktopWidget>
#include <QtWidgets/QMainWindow>
#include <iostream>
#include "cth_screeninfo.h"

using namespace std;

cth_ScreenInfo::cth_ScreenInfo(QWidget *parent, QScreen *screen) : QMainWindow(parent)
{
  curr_screen = screen;

  setWindowFlags(Qt::Widget | Qt::FramelessWindowHint); // invisible window
  setAttribute(Qt::WA_NoSystemBackground, true);
  setAttribute(Qt::WA_TranslucentBackground, true);
//  setAttribute(Qt::WA_PaintOnScreen);
}

cth_ScreenInfo::~cth_ScreenInfo()
{

   cout << size().width()  << "x" << 
           size().height() << "+" <<
           frameGeometry().x() << "+"  <<
           frameGeometry().y() << "+" << 
           curr_screen->geometry().width() << "+" <<
           curr_screen->geometry().height() << endl;
}

bool cth_ScreenInfo::event(QEvent *event)
{
   const bool ret_val = QMainWindow::event(event);
   if(!functionAfterShownCalled && event->type() == QEvent::Paint)
   {
      close();
   }
   return ret_val;
}


int main(int argc, char *argv[])
{
   int num_screens = 0;
   QApplication a(argc, argv);
   auto m = QApplication::desktop();

   foreach (QScreen *screen, QGuiApplication::screens()) 
   {
      auto desk = m->availableGeometry(num_screens);
      cth_ScreenInfo w(0,screen);
      w.move(desk.topLeft());
      w.showMaximized();
      ++num_screens;
      a.exec();
   }
}
