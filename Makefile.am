## Process this file with automake to produce Makefile.in
## In turn, Makefinle.in will be used by the configure script to generate the
## final Makefile.
SUBDIRS = src contrib

if SINGLE_POINT
SUBDIRS += data/offline/single_point/conf data/offline/single_point/params
DATAFILES_DIR = offline/single_point
BINFILE = inland-single_point
MODEL_INPUTDIR = input-single_point
else
SUBDIRS += data/offline/grid/conf data/offline/grid/params
DATAFILES_DIR = offline/grid
BINFILE = inland-grid
MODEL_INPUTDIR = input
endif

dist_doc_DATA = docs/README.grid docs/README.single_point docs/INSTALL.single_point docs/INSTALL.grid docs/guiarapido.txt docs/quickguide.txt
CLEANFILES = *.nc ibis_restart_control fort.* yearsrun.dat
DISTCLEANFILES = inland-single_point inland-grid conf params input output
EXTRA_DIST = include/inland_comage.h include/inland_compar.h include/inland_comsat.h include/inland_comtex.h

dev-symlinks: all
	@echo Creating symlinks to run the model from the current directory.
# first remove any existing links (useful when switching from 2d/1d model)
	for lnk in conf params input output ; do if [ -h $$lnk ]; then rm $$lnk ; fi ; done
# now create any missing links
	if [ ! -e $(BINFILE) ]; then ln -s "${srcdir}"/src/$(BINFILE) . ; fi
	if [ ! -e conf ] ; then ln -s "${srcdir}"/data/$(DATAFILES_DIR)/conf . ; fi
	if [ ! -e params ] ; then ln -s "${srcdir}"/data/$(DATAFILES_DIR)/params . ; fi
	if [ ! -e input ] ; then ln -s "../inland-data/"$(MODEL_INPUTDIR) input ; fi
	if [ ! -e "${srcdir}"/data/$(DATAFILES_DIR)/output ] ; then mkdir "${srcdir}"/data/$(DATAFILES_DIR)/output ; fi
	if [ ! -e output ] ; then ln -s "${srcdir}"/data/$(DATAFILES_DIR)/output output ; fi
	@echo Done.
