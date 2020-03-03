# The list of .edt, etc. files is contained in an Excel spreadsheet.
# This may change.  If you add a new one, make it the current one
#real file
current_xls=data_files_for_CTH_cluster_project_v37_all.xls
one_xls=data_files_for_CTH_cluster_project_v37_oneexp.xls
#symbolic link to it
curr_data=cth_data.xls

DEBUG_OR_NOT= -ggdb3
#DEBUG_OR_NOT= -O2
AUTOMAKE_OPTIONS= -Wno-portability

AM_CXXFLAGS = $(DEBUG_OR_NOT) -Wall
AM_CPPFLAGS = $(DEBUG_OR_NOT) -Wall -DDATADIR=\"$(pkgdatadir)/$(curr_data)\" -DVERSION=\"$(VERSION)\"
AM_CFLAGS =$(DEBUG_OR_NOT) -Wall -std=c99 -pg

LDADD = -lm -lboost_program_options
AM_LDFLAGS =

bin_PROGRAMS = cthgui cth_cluster cth_screeninfo
noinst_PROGRAMS = xlsrd test_xl 

dist_bin_SCRIPTS = cth_env run_cthgui run_cthterm cth_curve.R install_missing_R install_missing_octave cth_cells.pl cth_data.pl cth_xls2csv.pl merge_arch.m

dist_pkgdata_DATA = $(current_xls) $(one_xls) csplinedesx_$(VERSION).tar.gz glmfitpx_$(VERSION).tar.gz allexp_from_r.txt 100vago_8i_8lrm_8e.type 100vago_lrm_8_dend.type 100vago_i_8_dend.type 100vago_e_8_dend.type Archetypes_24_v1.type SwallowArchetypes_v1.type ALLswallowArchetypes_v1.type

doc_DATA = CTH_CLUSTERING.odt CTH_CLUSTERING.pdf ChangeLog COPYING LICENSE

BUILT_SOURCES = ui_cthgui.h qrc_cthgui.cpp moc_cthgui.cpp ui_helpbox.h moc_cthguidropedit.cpp moc_helpbox.cpp csplinedesx_$(VERSION).tar.gz glmfitpx_$(VERSION).tar.gz 

cthgui_code = main.cpp \
                 cthgui.cpp \
					  cth_impl.cpp \
					  cthguidropedit.cpp \
					  helpbox.cpp \
					  cthgui.ui \
					  cthgui.h \
					  cthguidropedit.h \
					  helpbox.h \
					  cthgui.qrc \
					  cthgui.ui \
					  helpbox.ui \
					  cthgui.pro \
                 cthgui.pro.shared \
					  cthicon.png 

cthgui_SOURCES = $(cthgui_code) $(BUILT_SOURCES)

CLEANFILES = ${BUILT_SOURCES}

cth_cluster_SOURCES = cth_cluster.cpp cth_cluster.h freexl.h freexl_internals.h freexl.c
cth_cluster_CXXFLAGS = -D_GLIBCXX_DEBUG $(DEBUG_OR_NOT)

xlsrd_SOURCES = xlsrd.cpp freexl.h freexl_internals.h freexl.c
xlsrd_CPPFLAGS = -D_REENTRANT -fPIC
xlsrd_CFLAGS = $(AM_CFLAGS)
test_xl_SOURCES = test_xl.c freexl.h freexl_internals.h freexl.c
text_xl_CFLAGS = -D_REENTRANT -fPIC $(AM_CFLAGS)

cthgui_CXXFLAGS = `pkg-config --cflags Qt5Gui Qt5Core Qt5Widgets` -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC ${DEFINES}
cthgui_LDFLAGS = `pkg-config --libs Qt5Gui Qt5Core Qt5Widgets`  -lGL -lpthread 

cth_screeninfo_SOURCES = cth_screeninfo.cpp cth_screeninfo.h moc_cth_screeninfo.cpp
cth_screeninfo_CXXFLAGS = `pkg-config --cflags Qt5Gui Qt5Core Qt5Widgets` -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC ${DEFINES}
cth_screeninfo_LDFLAGS = `pkg-config --libs Qt5Gui Qt5Core Qt5Widgets`  -lGL -lpthread 

moc_%.cpp: %.h
	@MOC@ -o $@ $<

ui_%.h: %.ui
	@UIC@ -o $@ $<

qrc_%.cpp: %.qrc
	@RCC@ -o $@ $<

#depends on everything in csplinedex dir
csplinedesx_$(VERSION).tar.gz: $(shell find csplinedesx -type f)
	sed -i -e 's/^Version:.*/Version: $(VERSION)/g' csplinedesx/DESCRIPTION
	R CMD build csplinedesx

#depends on everything in glmfitpx dir
glmfitpx_$(VERSION).tar.gz: $(shell find glmfitpx -type f)
	sed -i -e 's/^Version:.*/Version: $(VERSION)/g' glmfitpx/DESCRIPTION
	R CMD build glmfitpx

EXTRA_DIST = debian csplinedesx glmfitpx CTH_CLUSTERING.odt CTH_CLUSTERING.pdf

checkin_release:
	git add $(cthgui_code) $(cth_cluster_SOURCES) $(xlsrd_SOURCES) $(cth_screeninfo_SOURCES) $(current_xls) $(one_xls) $(EXTRA_DIST) Makefile.am configure.ac && git commit -m "Release cth files for version $(VERSION)"

checkpoint:
	git add $(cthgui_code) $(cth_cluster_SOURCES) $(xlsrd_SOURCES) $(cth_screeninfo_SOURCES) $(current_xls) $(one_xls) $(EXTRA_DIST) Makefile.am configure.ac && git commit --q -m "cth checkpoint files for version $(VERSION)"

checkpoint_withcomment:
	git add $(cthgui_code) $(cth_cluster_SOURCES) $(xlsrd_SOURCES) $(cth_screeninfo_SOURCES) $(current_xls) $(one_xls) $(EXTRA_DIST) Makefile.am configure.ac && git commit -q

# keep symbolic links to spreadsheet file current
install-data-hook:
	sed -i -e 's|_pkgdatadir_|$(pkgdatadir)|g' -e 's|_current_xls_|$(current_xls)|g' -e 's|_curr_data_|$(curr_data)|g' debian/cth-cluster.links


deb:
	@echo 'Making debian packages'
	make distdir &&\
	cd $(distdir) &&\
	dpkg-buildpackage -Zgzip -uc -us && \
	cd .. && \
	rm -rf $(distdir)

