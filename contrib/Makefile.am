contribdir = $(docdir)/contrib
dist_contrib_DATA = compare_0D-oct2010.ncl compare_0D-apr2012.ncl graphs_0D-apr2012.ncl open_0d.jnl convert_0d.jnl submit-inland.pbs

# if we are building grid model
# and --enable-highres_map was passed to configure, 
# we must build the inland-hrmap executable
# this requires that src/ is built first, which usually happens
if !SINGLE_POINT
if HIGHRES_MAP

CLEANFILES=*.mod *.oo
AM_CPPFLAGS = -I../src

bin_PROGRAMS = inland-write-hrmap
inland_write_hrmap_SOURCES = inland_write_hrmap.F90 ../src/inland_netcdfutils.F90
inland_write_hrmap.o: inland_write_hrmap.F90 #../src/inland_netcdfutils.o

endif
endif
