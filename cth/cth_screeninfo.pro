TEMPLATE = app
TARGET = cth_screeninfo

QT += core gui widgets

# Input
SOURCES += cth_screeninfo.cpp
HEADERS += cth_screeninfo.h

unix:!macx:!symbian: LIBS +=

lessThan(QT_MAJOR_VERSION, 5): error(This project requires Qt 5 or later)
MAKEFILE=Makefile_cthscreen.qt
