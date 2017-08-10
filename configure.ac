AC_INIT([inland], [1.0.0], [marcos.sanches@inpe.br])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_LANG([Fortran])
#AC_PROG_CC([gcc icc pgcc cc])

# Fortran-related macros
AC_PROG_FC([gfortran ifort pgf90 ftn],[Fortran 90])
AC_PROG_FC_C_O
AC_FC_SRCEXT(F90)
AC_FC_FREEFORM

# Special Fortran-90 macros
# You can find them in the autoconf macro archive at the address 
# http://www.gnu.org/software/autoconf-archive/. Their names are
# AX_F90_MODULE_EXTENSION and AX_F90_MODULE_FLAG.
#  Some tips presented at http://www.spinics.net/lists/ac/msg10297.html

# Find out how to add to the Fortran include path and store in ax_cv_f90_modflag
m4_include([m4/ax_f90_module_flag.m4])
AX_F90_MODULE_FLAG

# The module flag may be "-I " or so, and some make implementations
# swallow trailing space in variable settings, so we employ a little
# trick.
#AC_SUBST([empty], [""])
#AC_SUBST([FC_MODINC], [$ax_cv_f90_modflag'${empty}'])
ax_cv_f90_modflag="$(echo "$ax_cv_f90_modflag" | sed "s/ //g")"

# Define general macros
AC_DEFINE([STDOUT],[6],[Standard output unit])
AC_DEFINE([STDERR],[0],[Standard error unit])

# check if using cray's ftn compiler
AC_MSG_CHECKING([if we are using cray's ftn compiler])
if test "$FC" == "ftn"; then
  # ftn is cray's generic compiler for compute nodes 
  # ftn chooses compiler dependent on which PrgEnv-* module is loaded
  # set MYFC to specific compiler used based on PE_ENV - PATHSCALE and CRAY not supported yet
  if test "$PE_ENV" == "GNU"; then
    MYFC="gfortran"
  elif test "$PE_ENV" == "PGI"; then
    MYFC="pgf90"
  else
    MYFC="$PE_ENV (unknown)"
  fi
  AC_MSG_RESULT([yes, with $MYFC compiler])
else
  MYFC=$FC
  AC_MSG_RESULT([no])
fi

# define compiler-dependent preprocessor macros.
AC_MSG_CHECKING([which compiler are we using])
if test "$MYFC" == "gfortran"; then
  AC_DEFINE([GFORTRAN],[],[Define if using gcc/gfortran compiler.])
  FCFLAGS="-ffree-form -ffree-line-length-none ${FCFLAGS}"
  AC_MSG_RESULT([gcc/gfortran])
elif test "$MYFC" == "ifort"; then
  AC_DEFINE([IFORT],[],[Define if using intel/ifort compiler.])
  AC_DEFINE([IFORT_NOFLUSH_BUG],[],[Does ifort has the bug where advance=no does not flush once writing to STDOUT?])
  AC_MSG_RESULT([intel/ifort])
elif test "$MYFC" == "pgf90"; then
  AC_DEFINE([PGF90],[],[Define if using portland/pgf90 compiler.])
  AC_DEFINE([PGF90_NOFLUSH_BUG],[],[Does pgf90 has the bug where advance=no does not flush once writing to STDOUT?])
  AC_MSG_RESULT([portland/pgf90])
else
  AC_MSG_RESULT([unknown])
  #AC_MSG_WARNING([Unsupported FORTRAN compiler. Beware of unexpected compiler errors.])
fi

m4_include([m4/ac_fc_check_bounds.m4])
AC_FC_CHECK_BOUNDS()

# Enable using 'find' to locate library paths (slower)
AC_ARG_ENABLE([find_searches],
  [AS_HELP_STRING([--enable-find_searches],
     [Allow using the 'find' program to locate possible missing library paths. This is slower than using slocate, but is always up-to-date. @<:@disabled@:>@])],
  inland_lib_enable_find=true,inland_lib_enable_find=false)

AC_ARG_ENABLE([compiler_optimizations],
  [AS_HELP_STRING([--enable-compiler-optimizations],
     [Enable aggressive and often unreliable compiler optimizations that potentially can change the end results of a simulation. @<:@disabled@:>@])],
  inland_comp_optim_enable=true,inland_comp_optim_enable=false)

AC_ARG_ENABLE([intel_heap_arrays],
  [AS_HELP_STRING([--disable-intel-heap-arrays],
     [Disable intel compiler array allocation in heap. By default Intel compiler uses the stack to allocate arrays. In a global simulation, the minimum required stack size is 73728k (72M - ulimit -s 73728). If you are running a global simulation or if your stack size is too small and you disable this option, the model will die in a segfault. @<:@enabled@:>@])],
  inland_intel_heaparr_flag=false,inland_intel_heaparr_flag=true)

if ${inland_intel_heaparr_flag}; then
  AC_MSG_CHECKING([whether we need to add -heap-arrays compiler option])
  if test "${FC}" == "ifort"; then
    AC_MSG_RESULT([yes])
    AC_MSG_CHECKING([for maximum stack size])
    stacksz="$(ulimit -s 2>/dev/null)"
    if test "${stacksz}" == "$(echo "${stacksz}" | sed s/@<:@^0-9@:>@//g)"; then
      AC_MSG_RESULT([${stacksz} - $(echo "${stacksz}" | sed s/@<:@^0-9@:>@//g)])
    else
      AC_MSG_RESULT([absent])
      stacksz=""
    fi
    FCFLAGS="-heap-arrays $stacksz ${FCFLAGS}"
  else
    AC_MSG_RESULT([no])
  fi
fi
if ${inland_comp_optim_enable}; then
  AC_MSG_CHECKING([the compiler to which add optimization options])

  if test "$MYFC" == "gfortran"; then
    AC_MSG_RESULT([gfortran])
    # rip off default compiler flags
    FCFLAGS="$(echo "${FCFLAGS}" | sed "s/\(-g\|-O2\)//g")"
    FCFLAGS="-O3 -ffast-math -march=native -funroll-all-loops -fprefetch-loop-arrays -pipe -fPIC -m64 ${FCFLAGS}"
  elif test "$MYFC" == "pgf90"; then
    AC_MSG_RESULT([pgf90])
    # rip off default compiler flags
    FCFLAGS="$(echo "${FCFLAGS}" | sed "s/\(-g\|-O2\)//g")"
    FCFLAGS="-fast -Mipa=fast,inline -Mprefetch ${FCFLAGS}" #-Mfprelaxed -relax floating point precision, trade accuracy for speed.
  elif test "$MYFC" == "ifort"; then
    AC_MSG_RESULT([ifort])
    # rip off default compiler flags
    FCFLAGS="$(echo "${FCFLAGS}" | sed "s/\(-g\|-O2\)//g")"

    # on ifort, we can query the CPU to know the best instruction set to use
    AC_MSG_CHECKING([for /proc/cpuinfo])
    if test -e "/proc/cpuinfo"; then
      AC_MSG_RESULT([present])
      buff="$(grep "^flags" /proc/cpuinfo | head -n1)"
      axflag=""
      mtuneflag=""

      AC_MSG_CHECKING([for processor instruction set])
      # we don't know yet how would AVX instruction set show up on cpuinfo.
      if test ! -z "$(echo "${buff}" | grep " sse4_2\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSE4.2])
	axflag="-axSSE4.2"
	mtuneflag="-mtune=core2"
      elif test ! -z "$(echo "${buff}" | grep " sse4_1\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSE4.1])
	axflag="-axSSE4.1"
	mtuneflag="-mtune=core2"
      elif test ! -z "$(echo "${buff}" | grep " ssse3\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSSE3])
	axflag="-axSSSE3"
	mtuneflag="-mtune=core2"
      # must assure SSE3 instructions imply on core2 family processor!
      elif test ! -z "$(echo "${buff}" | grep " sse3\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSE3])
	axflag="-axSSE3"
	mtuneflag="-mtune=core2"
      elif test ! -z "$(echo "${buff}" | grep " sse2\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSE2])
	axflag="-axSSE2"
	if test ! -z "$(echo "${buff}" | grep " mmx\( \(\$\)")"; then
	  mtuneflag="-mtune=pentium4m"
	else
	  mtuneflag="-mtune=pentium4"
	fi
      elif test ! -z "$(echo "${buff}" | grep " sse\( \|\$\)")"; then
        AC_MSG_RESULT([up to SSE])
	mtuneflag="-mtune=pentiumpro"
      elif test ! -z "$(echo "${buff}" | grep " mmx\( \|\$\)")"; then
        AC_MSG_RESULT([up to MMX])
	mtuneflag="-mtune=pentium-mmx"
      fi
    else # no /proc/cpuinfo. can't probe CPU capabilities! blame macosx
      AC_MSG_RESULT([absent])
    fi

    FCFLAGS="${axflag} -ftz -ip -ipo ${mtuneflag} -O3 -xHost ${FCFLAGS}" #it's possible use -override-limits (biggest computer)
  else # don't know which compiler is this 
    AC_MSG_RESULT([unknown/unsupported])
  fi
fi

# openmp - tested for gfortran on linux+tupa, pgi on tupa
# should print warning when enable single_point
AC_OPENMP
if test "$OPENMP_FCFLAGS" != ""; then
     FCFLAGS="$OPENMP_FCFLAGS ${FCFLAGS}"
#   if test "$FC" == "ftn"; then
#     FC="ftn $OPENMP_FCFLAGS "
#   else
#     FCFLAGS="$OPENMP_FCFLAGS ${FCFLAGS}"
#  fi
fi


AC_ARG_WITH([netcdf],
  [AS_HELP_STRING([--with-netcdf],[Specify NetCDF custom location.])],
  [FCFLAGS="$FCFLAGS -I$withval/include"; LDFLAGS="$LDFLAGS -L$withval/lib"])

AC_ARG_WITH([hdf5],
  [AS_HELP_STRING([--with-hdf5],[Specify HDF5 custom location.])],
  [FCFLAGS="$FCFLAGS -I$withval/include"; LDFLAGS="$LDFLAGS -L$withval/lib"])

AC_ARG_WITH([curl],
  [AS_HELP_STRING([--with-curl],[Specify Curl library custom location.])],
  [FCFLAGS="$FCFLAGS -I$withval/include"; LDFLAGS="$LDFLAGS -L$withval/lib"])

# Different model modes (2D, 0D)
AC_ARG_ENABLE([single_point],
  [AS_HELP_STRING([--enable-single_point],
     [Compile the model as single-point (0D) version. @<:@disabled@:>@])],
  inland_single_point_compile=true,inland_single_point_compile=false)
# Define offline compile if this is the case
AM_CONDITIONAL([SINGLE_POINT],[$inland_single_point_compile])
if ${inland_single_point_compile}; then
 AC_DEFINE([SINGLE_POINT_MODEL],[],[Define if compiling single point (0D) model. Otherwise (for 2D/grid model) leave it undefined.])
fi

# agro model
AC_ARG_ENABLE([agro],
  [AS_HELP_STRING([--enable-agro],
     [Enable agro model. @<:@disabled@:>@])],
  inland_agro_compile=true,inland_agro_compile=false)

AM_CONDITIONAL([AGRO],[$inland_agro_compile])
if ${inland_agro_compile}; then
 AC_DEFINE([AGRO_MODEL],[],[Define if agro model is enabled.])
fi

# enable high-res map reading and writing (used in tiling scheme)
AC_ARG_ENABLE([highres_map],
  [AS_HELP_STRING([--enable-highres_map],
     [Enable high-res map reading and writing, used in tiling scheme. @<:@disabled@:>@])],
  inland_hrmap_compile=true,inland_hrmap_compile=false)

AM_CONDITIONAL([HIGHRES_MAP],[$inland_hrmap_compile])
if ${inland_hrmap_compile}; then
 AC_DEFINE([HRMAP],[],[Define if high-res map reading and writing is enabled.])
fi

# Check the NetCDF Library
function inland_netcdf_link_attempt() {
 if ! ${inland_single_point_compile}; then
  inland_netcdf_linkable=false
  if test ! -z "$FCFLAGS"; then
   AC_MSG_NOTICE([Using custom FCFLAGS: $FCFLAGS])
  fi
  if test ! -z "$LDFLAGS"; then
   AC_MSG_NOTICE([Using custom LDFLAGS: $LDFLAGS])
  fi
 
  for ncflibname in netcdf netcdff; do
   if ! $inland_netcdf_linkable; then
    unset ac_cv_lib_${ncflibname}___nf_inq_libvers
    AC_CHECK_LIB([${ncflibname}], [nf_inq_libvers],
                 AC_MSG_NOTICE([NetCDF needs no additional libraries!]);
                 inland_netcdf_linkable=true;
 		LIBS="$LIBS -l${ncflibname}")
   fi
   if ! $inland_netcdf_linkable; then
    unset ac_cv_lib_${ncflibname}__lhdf5_hl__lcurl___nf_inq_libvers
    AC_CHECK_LIB([${ncflibname} -lhdf5_hl -lcurl], [nf_inq_libvers],
            AC_MSG_NOTICE([NetCDF works with -lhdf5_hl -lcurl (HDF5 + libCurl)]);
 	   inland_netcdf_linkable=true;
 	   LIBS="$LIBS -l${ncflibname} -lhdf5_hl -lcurl"
   )
   fi
  done
else
 inland_netcdf_linkable=true
 AC_MSG_NOTICE([Single point (0D) model does not use NETCDF input/output.])
fi
}
 
# first attempt with provided LDFLAGS (-I / -L arguments)
inland_netcdf_link_attempt

# could not link. search the computer for a working netcdf lib
if ! $inland_netcdf_linkable; then 
 # function to get paths to possible netcdf libraries
 function inland_netcdf_getpath_attempt() {
  inland_libnetcdf_attempt="$1"
  if test ! -z "$inland_libnetcdf_attempt"; then
   inland_FCFLAGS_original=$FCFLAGS
   inland_LDFLAGS_original=$LDFLAGS
   for inland_ldpath in $inland_libnetcdf_attempt; do
    if ! $inland_netcdf_linkable; then
     inland_LIB_attempt=${inland_ldpath%/*}
     inland_INC_attempt=${inland_ldpath%/*/*}/include
     AC_MSG_NOTICE([trying $ax_cv_f90_modflag$inland_INC_attempt -L$inland_LIB_attempt])
     FCFLAGS="$ax_cv_f90_modflag$inland_INC_attempt $inland_FCFLAGS_original"
     LDFLAGS="-L$inland_LIB_attempt $inland_LDFLAGS_original"
     inland_netcdf_link_attempt
     if ! $inland_netcdf_linkable; then
      FCFLAGS="${inland_FCFLAGS_original}"
      LDFLAGS="${inland_LDFLAGS_original}"
     fi
    fi
   done
  fi
 }
 AC_MSG_NOTICE([NetCDF library not found. Trying to locate a working library.])
 AC_MSG_CHECKING([for slocate])
 inland_locate_path="$(which locate)"
 if test ! -z "$inland_locate_path"; then
  AC_MSG_RESULT([$inland_locate_path])
  inland_libnetcdf_candidates="$($inland_locate_path libnetcdf.a)"

  # try results in /usr, except /usr/src
  inland_netcdf_getpath_attempt "$(echo "$inland_libnetcdf_candidates" | grep ^/usr | grep -v ^/usr/src/)"

  if ! $inland_netcdf_linkable; then
   # try results in HOME directory and subdirectories.
   inland_netcdf_getpath_attempt "$(echo "$inland_libnetcdf_candidates" | grep ^$HOME/)"
  fi

  if ! $inland_netcdf_linkable; then
   # try results in /opt directory
   inland_netcdf_getpath_attempt "$(echo "$inland_libnetcdf_candidates" | grep ^/opt/)"
  fi

# if you want to remove the ones inside netcdf source tree:
# :: echo "$inland_libnetcdf_candidates" | grep ^/usr | grep -v "libsrc.*/\.libs/libnetcdf\.a\$"

  if ! $inland_netcdf_linkable; then
   # try the remaining results (get the first that works)
   inland_netcdf_getpath_attempt "$inland_libnetcdf_candidates"
  fi
 else
  AC_MSG_RESULT([no])
  if $inland_lib_enable_find; then
   AC_MSG_CHECKING(["for findutils' find"])
   inland_findutils_find_path="$(which find)"
   if test ! -z "$inland_findutils_find_path"; then
    AC_MSG_RESULT([$inland_findutils_find_path])

    AC_MSG_NOTICE([Searching for libnetcdf.a on /usr (except /usr/src)])
    # try results in /usr, except /usr/src
    inland_netcdf_getpath_attempt "$($inland_findutils_find_path /usr -name libnetcdf.a | grep -v ^/usr/src)"

    if ! $inland_netcdf_linkable; then
     AC_MSG_NOTICE([Searching for libnetcdf.a on /opt)])
     # try results in /opt directory
     inland_netcdf_getpath_attempt "$($inland_findutils_find_path /opt -name libnetcdf.a)"
    fi

    if ! $inland_netcdf_linkable; then
     AC_MSG_NOTICE([Searching for libnetcdf.a on $HOME)])
     # try results in HOME directory and subdirectories.
     inland_netcdf_getpath_attempt "$($inland_findutils_find_path $HOME -name libnetcdf.a)"
    fi

    if ! $inland_netcdf_linkable; then
     AC_MSG_NOTICE([Searching for libnetcdf.a on the whole disk (/))])
     # try the remaining results (get the first that works)
     inland_netcdf_getpath_attempt "$($inland_findutils_find_path / -name libnetcdf.a)"
    fi
   fi
   AC_MSG_NOTICE(["  Unable to link against NetCDF library and the slocate tool is not available 
to provide a quick means to search for possible locations of the netcdf library.
  Please install netcdf and provide appropriate LDFLAGS variable informing the correct paths for NetCDF's include and lib directory.
  You can try a more aggressive search if you have findutils' find and enable 
the 'find_searches' option in the commandline arguments to configure.
  For this, run this again with ./configure --enable-find_searches"])
  fi
 fi
fi

if ! $inland_netcdf_linkable; then 
 AC_MSG_ERROR([could not link against a NetCDF library.])
fi

# find out if we have netcdf-4 API, requires nc-config and netcdf >= 4 to test
# this might fail if found using the esoteric search above, instead of system path or --with-netcdf
 
dnl find nc-config location
unset ac_cv_path_NETCDF_NCCONFIG
if test "$with_netcdf" = "yes" -o "$with_netcdf" = "" ; then
  AC_PATH_PROG(NETCDF_NCCONFIG, nc-config, no)
else
  tmp_path="$with_netcdf/bin$PATH_SEPARATOR$with_netcdf"
  AC_PATH_PROG(NETCDF_NCCONFIG, nc-config, no, $tmp_path)
fi

dnl test nc-config
if test "$NETCDF_NCCONFIG" = "no" ; then
  echo "        did not find nc-config, some features may be missing"
  echo "        use --with-netcdf=/path/to/netcdf or add nc-config to PATH"
  NETCDF_NCCONFIG=
else

  dnl test version
  NETCDF_VERSION=`$NETCDF_NCCONFIG --version`
  AC_MSG_RESULT([ got netcdf version="$NETCDF_VERSION"])

  dnl test for netcdf-4
  AC_MSG_CHECKING([for netcdf-4 (and HDF5) support in libnetcdf])
      if test "x$($NETCDF_NCCONFIG --has-nc4)" = "xyes"; then :
        AC_MSG_RESULT([yes])
        AC_DEFINE([NETCDF4],[],[netcdf-4])
      else
        AC_MSG_RESULT([no])
      fi
fi


# Find out the file name extensions of Fortran module files.
m4_include([m4/ax_f90_module_extension.m4])
AX_F90_MODULE_EXTENSION
if test "x$ax_cv_f90_modext" != xunknown; then
  FC_MODEXT=$ax_cv_f90_modext
else
  FC_MODEXT=$ax_cv_f90_modext
fi
AC_SUBST([FC_MODEXT])

# Continuing normal ac file
AC_CONFIG_HEADERS([include/inland_config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
contrib/Makefile
data/offline/grid/conf/Makefile
data/offline/grid/params/Makefile
data/offline/single_point/conf/Makefile
data/offline/single_point/params/Makefile
])
AC_OUTPUT
