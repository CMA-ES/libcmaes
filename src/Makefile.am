
libcmaesdir = $(libdir)
libcmaesdatadir = $(datadir)
libcmaesincludedir = $(includedir)/libcmaes

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=libcmaes_config.h cmastrategy.cc esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc bipopcmastrategy.h bipopcmastrategy.cc covarianceupdate.h covarianceupdate.cc acovarianceupdate.h acovarianceupdate.cc pwq_bound_strategy.h pwq_bound_strategy.cc

include_HEADERS = cmastrategy.h esoptimizer.h esostrategy.h cmasolutions.h parameters.h cmaparameters.h cmastopcriteria.h ipopcmastrategy.h bipopcmastrategy.h covarianceupdate.h acovarianceupdate.h pwq_bound_strategy.h

AM_CPPFLAGS=-I$(EIGEN3_INC)
AM_CXXFLAGS=-Wall -g -fopenmp
