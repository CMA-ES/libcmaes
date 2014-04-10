
bin_PROGRAMS=simple_test test_functions
simple_test_SOURCES=simple-test.cc
test_functions_SOURCES=test-functions.cc

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC)
AM_CXXFLAGS=-g -Wall -std=c++11 -fopenmp
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
