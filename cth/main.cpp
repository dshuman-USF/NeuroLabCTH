#include <QtWidgets/QApplication>
#include "cthgui.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    CTHgui w;
    int cmd;
    int portnum = 6000;
    QString port;
    bool debug = false;

    const struct option opts[] = {
       {"d", no_argument, nullptr, 'd'},
       {"p", required_argument, nullptr, 'p'},
       { 0,0,0,0} };

   while ((cmd = getopt_long_only(argc, argv, "", opts, nullptr )) != -1)
   {
      switch (cmd)
      {
         case 'd':
            debug = true;
            break;

         case 'p':
            port = optarg;
            bool ok;
            if (port.length() == 0)
            {
               cout << "Port number is missing, aborting. . .\n";
               exit(1);
            }
            else
            {
               portnum = port.toInt(&ok,10);
               if (!ok)
               {
                  cout << "Port number is not a number, aborting. . .\n";
                  exit(1);
               }
            }
            break;
         default:
            break;
      }
   }
   w.setPort(portnum);
   w.setDebug(debug);
   w.show();
   return a.exec();
}
