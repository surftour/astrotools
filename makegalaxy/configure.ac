#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([MakeGalaxy], [1.0], [Patrik Jonsson <sunrise@familjenjonsson.org>])
AC_CONFIG_SRCDIR([prototypes.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_STDC

AC_CHECK_LIB([hdf5],[main],, AC_MSG_ERROR(HDF5 library not found))
AC_CHECK_HEADERS([hdf5.h],,AC_MSG_ERROR(HDF5 headers not found))
AC_CHECK_LIB([gslcblas],[main],, AC_MSG_ERROR(Gsl cblas library not found))
AC_CHECK_LIB([gsl],[main],, AC_MSG_ERROR(Gsl library not found))

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
