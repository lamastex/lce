
bin_PROGRAMS = tut
tut_SOURCES = tutorial.cpp \
				myApp.cpp myApp.hpp \
				myFrame.cpp myFrame.hpp \
				mscpp_pop_struct_listCtrl.cpp mscpp_pop_struct_listCtrl.hpp \
				mscpp_treebook.cpp mscpp_treebook.hpp \
				mscpp_migmat_grid.cpp mscpp_migmat_grid.hpp \
				mscpp_pop_dialogue.cpp mscpp_pop_dialogue.hpp

AM_CPPFLAGS = -I$(top_srcdir)/hudson_ms_lib
					-I$(top_builddir)/hudson_ms_lib

LIBMSADD = ../../hudson_ms_lib/.libs/libmscplusplus-0.0.a

tut_LDADD = $(LIBMSADD)
