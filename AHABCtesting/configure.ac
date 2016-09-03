
AC_PREREQ([2.63])
AC_INIT([ahabc], [0.1], [jenny.harlow@pg.canterbury.ac.nz])
AM_INIT_AUTOMAKE([no-define])

AC_CONFIG_SRCDIR([runs/make_run_pops.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

#check for gsl
# see http://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [],
	AC_MSG_ERROR([cannot find -lgslcblas])) # action if not found
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [],
	AC_MSG_ERROR([cannot find -lgsl])) # action if not found

AC_CHECK_LIB([sequence], [main], [],
    AC_MSG_ERROR([cannot find -lsequence])) # action if not found

AC_CHECK_LIB([cxsc], [main], [],
    AC_MSG_ERROR([cannot find -lcxsc]))

# This boost library was chosen as there was a bug-fix release of spirit here
# which is compatible with the code. Also, its from 2004.
# It would be nice to be able to use spirit v2, but that isn't available until
# boost 1.41 (released Nov 2009)
AX_BOOST_BASE([1.31])
#AX_BOOST_PROGRAM_OPTIONS
#BOOST_REQUIRE([1.31])
#BOOST_PROGRAM_OPTIONS

# Checks for header files.
# Checks for header files.
AC_CHECK_HEADERS([unistd.h])
AC_HEADER_STDBOOL
AC_CXX_HAVE_SSTREAM
AC_CXX_HAVE_STRSTREAM
AC_CXX_HAVE_LONG_LONG

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T


AC_LANG(C++)

AC_SEARCH_LIBS([m])
AC_CHECK_FUNCS([floor, pow, sqrt]) 

# Check and set compiler characteristics.

CFLAGS="$CFLAGS  -g -O0 -Wall -pedantic"
#CXXFLAGS="$CXXFLAGS -g -O0 -Wall -Wextra -Werror -pedantic"
CXXFLAGS="$CXXFLAGS -g -O0 -Wall -Wextra -Werror=implicit -Werror=missing-include-dirs -pedantic"


AC_CONFIG_FILES([Makefile
				runs/Makefile
				])
AC_OUTPUT
