AC_PREREQ([2.67])
AC_INIT([cth-cluster],[1.3.14],[dshuman@usf.edu])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR(cth_cluster.cpp)
AC_CONFIG_HEADERS([config.h])
: ${CPPFLAGS=""}
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_PROG_CXX
AC_HEADER_STDC
AC_PROG_LN_S

PKG_PROG_PKG_CONFIG

# Check for Qt libraries
PKG_CHECK_MODULES(QT, [Qt5Core, Qt5Gui, Qt5Widgets], [], [AC_MSG_ERROR([Qt libraries are required.])])

# Retrieve Qt compilation and linker flags
CPPFLAGS="`$PKG_CONFIG --cflags-only-I Qt5Core Qt5Gui Qt5Widgets` $CPPFLAGS"
LDFLAGS="`$PKG_CONFIG --libs-only-L Qt5Core Qt5Gui Qt5Widgets` $LDFLAGS"
LIBS="`$PKG_CONFIG --libs-only-l Qt5Core Qt5Gui Qt5Widgets`  $LIBS"

if ! `$PKG_CONFIG --atleast-version=5.7.0 Qt5Core`; then
	AC_MSG_ERROR([Qt 5.7.0 or greater is required.])
fi

AC_CHECK_PROGS(MOC, [moc])
AC_CHECK_PROGS(UIC, [uic])
AC_CHECK_PROGS(RCC, [rcc])

if test -z "$MOC" || test -z "$UIC" || test -z "$RCC"; then 
	AC_MSG_ERROR([Qt utility programs moc, uic, and rcc are required.])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
