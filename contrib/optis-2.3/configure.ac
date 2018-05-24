#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([optis], [2.3], [Vitor Fontes <vitor.fontes@ufv.br>])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([subdir-objects])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_SED
AC_PROG_LN_S

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset mkdir pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
