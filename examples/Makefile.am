bin_PROGRAMS=sample_code
sample_code_SOURCES=sample-code.cc

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC)
AM_CXXFLAGS=-g -Wall -std=c++11 -fopenmp
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
