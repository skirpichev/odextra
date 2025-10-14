# Process this file with autoconf to produce a configure script.

AC_INIT([odextra],[0.1],[skirpichev@gmail.com])
AC_CONFIG_SRCDIR(odextra.h)

AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE

AC_PROG_MAKE_SET

# Check for which system.
AC_CANONICAL_HOST

# Disable unnecessary libtool tests...
define([AC_LIBTOOL_LANG_CXX_CONFIG],[:])dnl
define([AC_LIBTOOL_LANG_F77_CONFIG],[:])dnl
define([AC_LIBTOOL_LANG_GCJ_CONFIG],[:])dnl

# Checks for programs.
AC_LANG(C)
AC_PROG_CC_C99
AC_PROG_LIBTOOL

# Checks for libraries...
AC_CHECK_LIB(gslcblas,main,,AC_MSG_ERROR([C-language BLAS is required.]))
AC_CHECK_LIB(gsl,main,,AC_MSG_ERROR([GNU Scientific library is required.]))

dnl Check whether we need -lm.
LT_LIB_M
LIBS="$LIBS $LIBM"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# From GSL configure.ac:
AH_BOTTOM([#define RETURN_IF_NULL(x) if (!x) { return ; }
])
AH_BOTTOM([/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif])

AC_CONFIG_FILES([makefile doc/makefile tests/makefile])
AC_OUTPUT
