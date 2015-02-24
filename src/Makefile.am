
libcmaesdir = $(libdir)
libcmaesdatadir = $(datadir)
libcmaesincludedir = $(includedir)/libcmaes

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=libcmaes_config.h cmaes.h eo_matrix.h cmastrategy.cc esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc bipopcmastrategy.h bipopcmastrategy.cc covarianceupdate.h covarianceupdate.cc acovarianceupdate.h acovarianceupdate.cc vdcmaupdate.h vdcmaupdate.cc pwq_bound_strategy.h pwq_bound_strategy.cc eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h pli.h errstats.cc errstats.h contour.h

nobase_libcmaesinclude_HEADERS = cmaes.h opti_err.h eo_matrix.h cmastrategy.h esoptimizer.h esostrategy.h cmasolutions.h parameters.h cmaparameters.h cmastopcriteria.h ipopcmastrategy.h bipopcmastrategy.h covarianceupdate.h acovarianceupdate.h vdcmaupdate.h pwq_bound_strategy.h eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h errstats.h pli.h contour.h

if HAVE_SURROG
libcmaes_la_SOURCES += surrcmaes.h surrogatestrategy.cc surrogatestrategy.h surrogates/rankingsvm.hpp surrogates/rsvm_surr_strategy.hpp
nobase_libcmaesinclude_HEADERS += surrcmaes.h surrogatestrategy.h surrogates/rankingsvm.hpp surrogates/rsvm_surr_strategy.hpp
endif

AM_CPPFLAGS=-I$(EIGEN3_INC)
AM_CXXFLAGS=-Wall -Wextra -g
if !HAVE_CLANG
AM_CXXFLAGS += -fopenmp
endif

uninstall-hook:
	if test -d ${libcmaesincludedir}; then rm -rf ${libcmaesincludedir}; fi
