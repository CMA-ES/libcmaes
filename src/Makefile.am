
libcmaesdir =$(libdir)

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=cmaregular.h cmastrategy.h cmastrategy.cc covarianceupdate.h esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc

AM_CPPFLAGS=-I$(EIGEN3_INC)
AM_CXXFLAGS=-Wall -g -std=c++11 -fopenmp
