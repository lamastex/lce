#append flags
#ACLOCAL_AMFLAGS = -I m4 --install ${ACLOCAL_FLAGS}


#flag for config.h include
CONFIG_F=-I$(top_srcdir) -I$(top_builddir)



bin_PROGRAMS = 	testParameterParticleSet \
				testParameterParticleAnalysis \
				testSummaryStatistic \
				testAHABCValues \
				testArgConverter \
				testArgFactory \
				testMicrosatLocus \
				testMultiLociMicrosat \
				testPartitionedArgSeqGen \
				testMixtureTest1D \
				testMixtureTest2D \
				testPolyDataSequenceMaker \
				testPolysitesBasic \
				testPolysitesSmall \
				testSimdataBasic \
				testMicrosatBasic \
				testMultiLociMicrosatPopStats \
				testMicrosatSetPopStats \
				testPolysitesTiming \
				testSimdataTiming \
				testMultiLociPolyTablePopStats \
				testLinReg \
				testReadNucData
				
testParameterParticleSet_SOURCES = 		testParameterParticleSet.cpp

testParameterParticleAnalysis_SOURCES = testParameterParticleAnalysis.cpp

testSummaryStatistic_SOURCES = 			testSummaryStatistic.cpp

testAHABCValues_SOURCES = 				testAHABCValues.cpp

testArgConverter_SOURCES = 				testArgConverter.cpp \
										make_test_pops.cpp make_test_pops.hpp

testArgFactory_SOURCES = 				testArgFactory.cpp \
										make_test_pops.cpp make_test_pops.hpp

testMicrosatLocus_SOURCES =		 		testMicrosatLocus.cpp

testMultiLociMicrosat_SOURCES = 		testMultiLociMicrosat.cpp

testPartitionedArgSeqGen_SOURCES = 		testPartitionedArgSequenceGenerator.cpp \
										make_test_pops.cpp make_test_pops.hpp

testMixtureTest1D_SOURCES = 			testMixtureTest1D.cpp

testMixtureTest2D_SOURCES = 			testMixtureTest2D.cpp

testPolysitesBasic_SOURCES = 			testPolysitesBasic.cpp

testPolyDataSequenceMaker_SOURCES = 	testPolyDataSequenceMaker.cpp

testPolysitesSmall_SOURCES = 			testPolysitesSmall.cpp

testSimdataBasic_SOURCES = 				testSimdataBasic.cpp

testMicrosatBasic_SOURCES = 			testMicrosatBasic.cpp \
										make_test_pops.cpp make_test_pops.hpp

testMultiLociMicrosatPopStats_SOURCES = 	testMultiLociMicrosatPopStats.cpp \
										make_test_pops.cpp make_test_pops.hpp

testMicrosatSetPopStats_SOURCES = 			testMicrosatSetPopStats.cpp \
										make_test_pops.cpp make_test_pops.hpp

testPolysitesTiming_SOURCES = 			testPolysitesTiming.cpp

testSimdataTiming_SOURCES = 			testSimdataTiming.cpp

testMultiLociPolyTablePopStats_SOURCES =	testMultiLociPolyTablePopStats.cpp \
										make_test_pops.cpp make_test_pops.hpp

testLinReg_SOURCES = 					testLinReg.cpp

testReadNucData_SOURCES = 				testReadNucData.cpp

EXTRA_DIST = modelsites2.txt \
		tree.newick

#basic am compiler flags
AM_CPPFLAGS = $(CONFIG_F) @SRC_SUBS_INCLUDES@ $(BOOST_CPPFLAGS)

#all the flags and make sure we have the AM_CPPFLAGS too or it will be overridden	
#test_CPPFLAGS =  $(AM_CPPFLAGS)

AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 

#additional library for every program
LDADD = ${top_builddir}/src/libmct.a

#program-specific additions to LDADD - this overrides LDADD so I make that part of the prog-specific thing
#important - add libmct.a (in LDADD) first, else libmct.a does not link properly to libmscplusplus
testArgConverter_LDADD = ${LDADD} -lmscplusplus-0.0 

testArgFactory_LDADD = ${LDADD} -lmscplusplus-0.0 

testMultiLociMicrosat_LDADD = ${LDADD} -lmscplusplus-0.0 

testPartitionedArgSeqGen_LDADD = ${LDADD} -lmscplusplus-0.0 

testPolysitesBasic_LDADD = ${LDADD} -lmscplusplus-0.0 

testPolysitesSmall_LDADD = ${LDADD} -lmscplusplus-0.0 

testSimdataBasic_LDADD = ${LDADD} -lmscplusplus-0.0 

testMicrosatBasic_LDADD = ${LDADD} -lmscplusplus-0.0 

testMicrosatSetPopStats_LDADD = ${LDADD} -lmscplusplus-0.0 

testMultiLociMicrosatPopStats_LDADD = ${LDADD} -lmscplusplus-0.0 

testPolysitesTiming_LDADD = ${LDADD} -lmscplusplus-0.0 

testSimdataTiming_LDADD = ${LDADD} -lmscplusplus-0.0 

testMultiLociPolyTablePopStats_LDADD = ${LDADD} -lmscplusplus-0.0 

testLinReg_LDADD = ${LDADD} -lmscplusplus-0.0 

testReadNucData_LDADD = ${LDADD} -lmscplusplus-0.0 

LIBS += $(BOOST_PROGRAM_OPTIONS_LIBS)
