
libcmaesdir = $(libdir)
libcmaesdatadir = $(datadir)
libcmaesincludedir = $(includedir)/libcmaes

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=libcmaes_config.h cmaes.h eo_matrix.h cmastrategy.cc esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc bipopcmastrategy.h bipopcmastrategy.cc covarianceupdate.h covarianceupdate.cc acovarianceupdate.h acovarianceupdate.cc pwq_bound_strategy.h pwq_bound_strategy.cc eigenmvn.h candidate.h genopheno.h noboundstrategy.h

libcmaesinclude_HEADERS = cmaes.h eo_matrix.h cmastrategy.h esoptimizer.h esostrategy.h cmasolutions.h parameters.h cmaparameters.h cmastopcriteria.h ipopcmastrategy.h bipopcmastrategy.h covarianceupdate.h acovarianceupdate.h pwq_bound_strategy.h eigenmvn.h candidate.h genopheno.h noboundstrategy.h

AM_CPPFLAGS=-I$(EIGEN3_INC)
AM_CXXFLAGS=-Wall -g -fopenmp
