PREFIX=$(HOME)#
#CXXINC=-I$(PREFIX)/include -L/usr/local/lib -I.#
CXXINC=-I$(PREFIX)/include -I/usr/local/include -L$(PREFIX)/lib -L/usr/local/lib -I.#
#CXXINC=-I$(PREFIX)/include -L$(PREFIX)/lib -I.#
CXX=c++
CC=cc
DEBUG=-DNDEBUG
OPT=-O2
CXXFLAGS=-Wall -W -ansi -pedantic -fabi-version=0 $(OPT) $(DEBUG) $(CXXINC)
CFLAGS=-Wall -W -ansi -pedantic  $(OPT) $(DEBUG) $(CXXINC)
LIBS=-lsequence -lgsl -lgslcblas #-lm

TARGET_PARAMEST=ParamEst
TARGET_TAJDPVALS1=TajDPvals1
TARGET_TAJDPVALSM=TajDPvalsM
ALL_TARGETS= TARGET_PARAMEST TARGET_TAJDPVALS1 TARGET_TAJDPVALSM

SUPP_OBJECTS=ebc_moves.o ebc_summaries.o ebc_output.o ebc_sfsinteg.o ebc_sfsparticles.o ebc_params.o ebc_process_options.o ebc_sfsdata.o ebc_fseq.o ebc_fseq_general_graph.o ebc_fxzseq_named_graph.o ebc_fxzseq_efficient_graph.o ebc_fseq_efficient_graph.o ebc_precomp_TsPs.o getopt.o getopt1.o ebc_sfstypes.o

OBJECTS_PARAMEST=$(SUPP_OBJECTS) $(TARGET_PARAMEST).o

OBJECTS_TAJDPVALS1=$(SUPP_OBJECTS) $(TARGET_TAJDPVALS1).o

OBJECTS_TAJDPVALSM=$(SUPP_OBJECTS) $(TARGET_TAJDPVALSM).o

ALL_OBJECTS=$(SUPP_OBJECTS) $(TARGET_PARAMEST).o $(TARGET_TAJDPVALS1).o $(TARGET_TAJDPVALSM).o

all: $(ALL_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET_PARAMEST) $(OBJECTS_PARAMEST) $(LIBS)
	@echo made $(TARGET_PARAMEST)
	$(CXX) $(CXXFLAGS) -o $(TARGET_TAJDPVALS1) $(OBJECTS_TAJDPVALS1) $(LIBS)
	@echo made $(TARGET_TAJDPVALS1)
	$(CXX) $(CXXFLAGS) -o $(TARGET_TAJDPVALSM) $(OBJECTS_TAJDPVALSM) $(LIBS)
	@echo made $(TARGET_TAJDPVALSM)


paramest: $(OBJECTS_PARAMEST)
	$(CXX) $(CXXFLAGS) -o $(TARGET_PARAMEST) $(OBJECTS_PARAMEST) $(LIBS)
	@echo made $(TARGET_PARAMEST)

tajdpvals1: $(OBJECTS_TAJDPVALS1)
	$(CXX) $(CXXFLAGS) -o $(TARGET_TAJDPVALS1) $(OBJECTS_TAJDPVALS1) $(LIBS)
	@echo made $(TARGET_TAJDPVALS1)

tajdpvalsM: $(OBJECTS_TAJDPVALSM)
	$(CXX) $(CXXFLAGS) -o $(TARGET_TAJDPVALSM) $(OBJECTS_TAJDPVALSM) $(LIBS)
	@echo made $(TARGET_TAJDPVALSM)

clean:
	rm -f $(ALL_TARGETS) $(ALL_OBJECTS)

dist:
	mkdir ebc_MT_project
	mkdir ebc_MT_project/work
	cp *.cpp *.hpp *.c *.h Makefile ebc_MT_project
	cp work/*moves ebc_MT_project/work/
	cp work/run.sh ebc_MT_project/work/
	cp work/sfsms.1 ebc_MT_project/work/
	tar zcvf ebc_MT_project.tgz ebc_MT_project/
	rm -rf ebc_MT_project

ALL_HEADERS=ebc_fseq.hpp ebc_fseq_general_graph.hpp ebc_fxzseq_named_graph.hpp ebc_fxzseq_efficient_graph.hpp ebc_fseq_efficient_graph.hpp ebc_moves.hpp ebc_output.hpp ebc_params.hpp  ebc_process_options.hpp ebc_sfsdata.hpp ebc_sfsinteg.hpp ebc_sfsparticles.hpp	ebc_summaries.hpp ebc_precomp_TsPs.hpp ebc_sfstypes.hpp ebc_graphtypes.hpp

TajDPvals1.o: TajDPvals1.cpp $(ALL_HEADERS)
TajDPvalsM.o: TajDPvalsM.cpp $(ALL_HEADERS)
ParamEst.o: ParamEst.cpp $(ALL_HEADERS)
ebc_fseq.o: ebc_sfstypes.hpp
ebc_fseq_general_graph.o: ebc_sfstypes.hpp ebc_graphtypes.hpp
ebc_fxzseq_named_graph.o: ebc_fseq_general_graph.hpp ebc_sfstypes.hpp ebc_graphtypes.hpp
ebc_fxzseq_efficient_graph.o ebc_fseq_efficient_graph.o: ebc_fseq_general_graph.hpp ebc_fxzseq_named_graph.hpp ebc_sfstypes.hpp ebc_graphtypes.hpp
ebc_moves.o: ebc_moves.cpp ebc_moves.hpp
ebc_output.o: ebc_output.cpp ebc_output.hpp
ebc_params.o: ebc_params.cpp ebc_params.hpp
ebc_process_options.o: ebc_process_options.cpp ebc_process_options.hpp getoptFix.h
ebc_sfsdata.o: ebc_sfsdata.cpp ebc_sfsdata.hpp
ebc_sfsinteg.o: ebc_sfsinteg.cpp ebc_sfsinteg.hpp
ebc_sfsparticles.o: ebc_sfsparticles.cpp ebc_sfsparticles.hpp
ebc_summaries.o: ebc_summaries.cpp ebc_summaries.hpp
ebc_sfstypes.o: ebc_sfstypes.cpp ebc_sfstypes.hpp
ebc_precomp_TsPs.o: ebc_precomp_TsPs.cpp ebc_precomp_TsPs.hpp

#down here, one can specify particular file dependencies,
#for example, if main.o requires models.hpp and sfsl.hpp,
#the line below will force main.o to be recompiled if either
#header changes

#main.o: main.cpp main.hpp ebc_sfstypes.hpp ebc_output.hpp ebc_sfsinteg.hpp ebc_moves.hpp ebc_summaries.hpp ebc_params.hpp ebc_process_options.hpp ebc_sfsdata.hpp ebc_fseq.hpp

#ebc_sfstypes.o: ebc_sfstypes.hpp

