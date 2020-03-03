#-------------------------------------------------
#
# Project created by QtCreator 2014-06-27T11:41:03
#
#-------------------------------------------------

QT += core gui widgets

TARGET = cthgui
TEMPLATE = app


SOURCES += main.cpp\
           cthgui.cpp \
           cth_impl.cpp \
    helpbox.cpp \
    cthguidropedit.cpp

HEADERS += cthgui.h \
    helpbox.h \
    cthguidropedit.h

FORMS += cthgui.ui \
    helpbox.ui

unix:!macx:!symbian: LIBS +=

OTHER_FILES +=

RESOURCES += cthgui.qrc

DEFINES += VERSION=\\\"1.3.9\\\"


# If you run qtcreator, it will clobber the autotools
# Makefile.  This causes qmake to output a Makefile
# with this name for development within the 
# qtcreator gui.
MAKEFILE=Makefile.qt
