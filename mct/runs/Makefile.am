#append flags
#ACLOCAL_AMFLAGS = -I m4 --install ${ACLOCAL_FLAGS}


#flag for config.h include
CONFIG_F=-I$(top_srcdir) -I$(top_builddir)



bin_PROGRAMS = 	runAABC_Simdata \
				runAABC_Simdata_FixObserved \
				runSample_Simdata \
				runSample_Polysites \
				runSample_Microsat \
				runSampleWithJiggle_Microsat \
				runSample_AHABC_nucs \
				runSample_AHABC_microsats
				
				
runAABC_Simdata_SOURCES = 					runAABC_Simdata.cpp

runAABC_Simdata_FixObserved_SOURCES =		runAABC_Simdata_FixObserved.cpp

runSample_Simdata_SOURCES =					runSample_Simdata.cpp

runSample_Polysites_SOURCES =				runSample_Polysites.cpp

runSample_Microsat_SOURCES =				runSample_Microsat.cpp \
											make_run_pops.cpp make_run_pops.hpp

runSampleWithJiggle_Microsat_SOURCES =		runSampleWithJiggle_Microsat.cpp \
											make_run_pops.cpp make_run_pops.hpp

runSample_AHABC_nucs_SOURCES =				runSample_AHABC_nucs.cpp \
											make_run_pops.cpp make_run_pops.hpp

runSample_AHABC_microsats_SOURCES =			runSample_AHABC_microsats.cpp \
											make_run_pops.cpp make_run_pops.hpp

EXTRA_DIST =

#basic am compiler flags
AM_CPPFLAGS = $(CONFIG_F) @SRC_SUBS_INCLUDES@ $(BOOST_CPPFLAGS)

#all the flags and make sure we have the AM_CPPFLAGS too or it will be overridden	
#test_CPPFLAGS =  $(AM_CPPFLAGS)

AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 

#additional library for every program
LDADD = ${top_builddir}/src/libmct.a

#program-specific additions to LDADD - this overrides LDADD so I make that part of the prog-specific thing
#important - add libmct.a (in LDADD) first
runAABC_Simdata_LDADD = ${LDADD} -lmscplusplus-0.0 

runAABC_Simdata_FixObserved_LDADD = ${LDADD} -lmscplusplus-0.0 

runSample_Simdata_LDADD = ${LDADD} -lmscplusplus-0.0 

runSample_Polysites_LDADD = ${LDADD} -lmscplusplus-0.0 

runSample_Microsat_LDADD = ${LDADD} -lmscplusplus-0.0 

runSampleWithJiggle_Microsat_LDADD = ${LDADD} -lmscplusplus-0.0 

runSample_AHABC_nucs_LDADD = ${LDADD} -lmscplusplus-0.0 

runSample_AHABC_microsats_LDADD = ${LDADD} -lmscplusplus-0.0 

LIBS += $(BOOST_PROGRAM_OPTIONS_LIBS)
