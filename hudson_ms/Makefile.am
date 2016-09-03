#conditional compilation - see configure.ac
# only try to make gui if we have wxWidgets library
if COND_WX
   MAYBE_GUI = gui
endif
     
SUBDIRS = hudson_ms_lib test $(MAYBE_GUI)
EXTRA_DIST = bootstrap
ACLOCAL_AMFLAGS = -I m4 -install
#ACLOCAL_AMFLAGS = -I m4 -install ${ACLOCAL_FLAGS}
#AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)
