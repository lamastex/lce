
AC_PREREQ([2.63])

AC_INIT([hudson_ms], [0.0], [r.sainudiin@math.canterbury.ac.nz])

AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE[-Wall -Werror]

AC_CONFIG_SRCDIR([test/make_test_pops.cpp])


AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# C++ language
#AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
#AC_PROG_RANLIB


# Checks for libraries.

# Use the wx-supplied m4 macros
# see http://wiki.wxwidgets.org/Autoconf
AM_OPTIONS_WXCONFIG
reqwx=2.9.0
AM_PATH_WXCONFIG($reqwx, wxWin=1)
if test "$wxWin" != 1; then
	AC_MSG_WARN([
		wxWidgets must be installed on your system to compile the gui.
 
		To compile the gui, please check that wx-config is in path, the directory
		where wxWidgets libraries are installed (returned by
		'wx-config --libs' or 'wx-config --static --libs' command)
		is in LD_LIBRARY_PATH or equivalent variable and
		wxWidgets version is $reqwx or above.
		])
	
fi
 
CPPFLAGS="$CPPFLAGS $WX_CPPFLAGS"
CXXFLAGS="$CXXFLAGS $WX_CXXFLAGS_ONLY"
CFLAGS="$CFLAGS $WX_CFLAGS_ONLY"
LIBS="$LIBS $WX_LIBS"


#check for gsl
# see http://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

#AC_CHECK_LIB([sequence], [main], [],
 #   AC_MSG_ERROR([cannot find -lsequence]))

# This boost library was chosen as there was a bug-fix release of spirit here
# which is compatible with the code. Also, its from 2004.
# It would be nice to be able to use spirit v2, but that isn't available until
# boost 1.41 (released Nov 2009)
AX_BOOST_BASE([1.31])
#AX_BOOST_PROGRAM_OPTIONS
#BOOST_REQUIRE([1.31])
#BOOST_PROGRAM_OPTIONS


# Checks for header files.

# Checks for typedefs and structures
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Check and set compiler characteristics.
#CC_CHECK_CFLAGS_APPEND([-Wall -Werror])

# wxWidgets fails to compile with -pedantic unless you ignore long longs
# see http://trac.wxwidgets.org/ticket/8517
CFLAGS="$CFLAGS  -g -O0 -Wall -pedantic"
CXXFLAGS="$CXXFLAGS -g -O0 -Wall -pedantic -Wno-long-long"

# Checks for library functions.



# Define these substitions here to keep all version information in one place.
# For information on how to properly maintain the library version information,
# refer to the libtool manual, section "Updating library version information":
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
#AC_SUBST([MSCPLUSPLUS_SO_VERSION], [0:0:0])
#AC_SUBST([MSCPLUSPLUS_API_VERSION], [0.0])


#conditional automake, only try to make gui if we have wxWidgets library
AM_CONDITIONAL([COND_WX], [test "$wxWin" == 1])
     
AC_CONFIG_FILES([Makefile
				test/Makefile 
				gui/Makefile
				gui/src/Makefile
				])

# the library has its own configure that outputs its own Makefile
AC_CONFIG_SUBDIRS([hudson_ms_lib])
AC_OUTPUT
