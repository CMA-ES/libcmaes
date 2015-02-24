AUTOMAKE_OPTIONS = subdir-objects
bin_PROGRAMS=simple_test edm
if HAVE_GFLAGS
bin_PROGRAMS += test_functions tpa_tests
if HAVE_SURROG
bin_PROGRAMS += surr_tests
surr_tests_SOURCES=surr-tests.cc
endif
test_functions_SOURCES=test-functions.cc
tpa_tests_SOURCES=tpa-tests.cc
endif
simple_test_SOURCES=simple-test.cc
edm_SOURCES=edm.cc
if HAVE_BBOB
libbbobdir = $(libdir)
libbbob_LIBRARIES=libbbob.a
libbbob_a_SOURCES= bbob.v13.09/c/benchmarksdeclare.c bbob.v13.09/c/benchmarkshelper.c bbob.v13.09/c/benchmarks.c bbob.v13.09/c/benchmarksnoisy.c bbob.v13.09/c/fgeneric.c bbob.v13.09/c/dirOK.c 
bin_PROGRAMS += bbobexperiment
bbobexperiment_SOURCES=bbobexperiment.cc
endif

if HAVE_GTEST
TESTS = $(check_PROGRAMS)
check_PROGRAMS = ut_pwqbounds ut_errstats ut_scaling
ut_pwqbounds_SOURCES=ut-pwqbounds.cc
ut_errstats_SOURCES=ut-errstats.cc
ut_scaling_SOURCES=ut-scaling.cc
endif

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC) $(GFLAGS_CFLAGS)
AM_CXXFLAGS=-g -Wall
if !HAVE_CLANG
AM_CXXFLAGS += -fopenmp
endif
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
if HAVE_BBOB
AM_CPPFLAGS += -Ibbob.v13.09/c
AM_LDFLAGS += -L.
LDADD += -lbbob
endif
if HAVE_GTEST
AM_LDFLAGS += -lgtest -lgtest_main
endif

if HAVE_BBOB
clean-local:
	rm -rf bbobc.tar.gz bbob.v13* *_bbob .deps/benchmarks* .deps/dir* .deps/fgeneric* .deps/bbob*
endif
