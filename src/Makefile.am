
libcmaesdir = $(libdir)
libcmaesdatadir = $(datadir)
libcmaesincludedir = $(includedir)/libcmaes

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=libcmaes_config.h cmaes.h eo_matrix.h cmastrategy.cc esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc bipopcmastrategy.h bipopcmastrategy.cc covarianceupdate.h covarianceupdate.cc acovarianceupdate.h acovarianceupdate.cc vdcmaupdate.h vdcmaupdate.cc pwq_bound_strategy.h pwq_bound_strategy.cc eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h surrogatestrategy.cc surrogatestrategy.h surrogates/rankingsvm.hpp surrogates/surr_rsvm_strategy.hpp pli.h errstats.cc errstats.h contour.h

libcmaesinclude_HEADERS = cmaes.h eo_matrix.h cmastrategy.h esoptimizer.h esostrategy.h cmasolutions.h parameters.h cmaparameters.h cmastopcriteria.h ipopcmastrategy.h bipopcmastrategy.h covarianceupdate.h acovarianceupdate.h vdcmaupdate.h pwq_bound_strategy.h eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h surrogatestrategy.h pli.h contour.h

AM_CPPFLAGS=-I$(EIGEN3_INC)
AM_CXXFLAGS=-Wall -Wextra -g -fopenmp

uninstall-hook:
	if test -d ${libcmaesincludedir}; then rm -rf ${libcmaesincludedir}; fi
