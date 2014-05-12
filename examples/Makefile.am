bin_PROGRAMS=sample_code sample_code_genopheno sample_code_pfunc
sample_code_SOURCES=sample-code.cc
sample_code_genopheno_SOURCES=sample-code-genopheno.cc
sample_code_pfunc_SOURCES=sample-code-pfunc.cc

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC)
AM_CXXFLAGS=-g -Wall -fopenmp
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
