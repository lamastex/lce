bin_PROGRAMS = test_pops test_rand test_segsites test_params test_args
test_pops_SOURCES = test_pops.cpp test_pops.hpp \
					test_pops_subs.cpp test_pops_subs.hpp \
					make_test_pops.cpp make_test_pops.hpp
test_rand_SOURCES = test_rand.cpp test_rand.hpp
test_segsites_SOURCES = test_segsites.cpp test_segsites.hpp \
						make_test_pops.cpp make_test_pops.hpp
test_params_SOURCES = test_params.cpp test_params.hpp \
						make_test_pops.cpp make_test_pops.hpp
test_args_SOURCES = test_args.cpp test_args.hpp \
						make_test_pops.cpp make_test_pops.hpp
						
AM_CPPFLAGS = -I$(top_srcdir)/hudson_ms_lib
					-I$(top_builddir)/hudson_ms_lib
#example_LDFLAGS = -static

LIBMSADD = ../hudson_ms_lib/.libs/libmscplusplus-0.0.a



test_pops_LDADD = $(LIBMSADD)
test_rand_LDADD = $(LIBMSADD)
test_segsites_LDADD = $(LIBMSADD)
test_params_LDADD = $(LIBMSADD)
test_args_LDADD = $(LIBMSADD)
