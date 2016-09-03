#append flags
#ACLOCAL_AMFLAGS = -I m4 --install ${ACLOCAL_FLAGS}


#flag for config.h include
CONFIG_F=-I$(top_srcdir) -I$(top_builddir)



bin_PROGRAMS = 	runNucScenarios
				
				

runNucScenarios_SOURCES =				runNucScenarios.cpp \
										make_run_pops.cpp make_run_pops.hpp

EXTRA_DIST =

#basic am compiler flags
AM_CPPFLAGS = $(CONFIG_F) $(BOOST_CPPFLAGS)

#all the flags and make sure we have the AM_CPPFLAGS too or it will be overridden	
#test_CPPFLAGS =  $(AM_CPPFLAGS)

AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 

#additional library for every program
LDADD = -lmct

#program-specific additions to LDADD - this overrides LDADD so I make that part of the prog-specific thing
#important - add libmct.a (in LDADD) first
runNucScenarios_LDADD = ${LDADD} -lmscplusplus-0.0 


LIBS += $(BOOST_PROGRAM_OPTIONS_LIBS)
