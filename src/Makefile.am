
libcmaesdir = $(libdir)
libcmaesdatadir = $(datadir)
libcmaesincludedir = $(includedir)

libcmaes_LTLIBRARIES=libcmaes.la
libcmaes_la_SOURCES=libcmaes_config.h cmaes.h eo_matrix.h cmastrategy.cc esoptimizer.h esostrategy.h esostrategy.cc cmasolutions.h cmasolutions.cc parameters.h cmaparameters.h cmaparameters.cc cmastopcriteria.h cmastopcriteria.cc ipopcmastrategy.h ipopcmastrategy.cc bipopcmastrategy.h bipopcmastrategy.cc covarianceupdate.h covarianceupdate.cc acovarianceupdate.h acovarianceupdate.cc vdcmaupdate.h vdcmaupdate.cc pwq_bound_strategy.h pwq_bound_strategy.cc eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h pli.h errstats.cc errstats.h contour.h

nobase_libcmaesinclude_HEADERS = ../include/libcmaes/cmaes.h ../include/libcmaes/opti_err.h ../include/libcmaes/eo_matrix.h ../include/libcmaes/cmastrategy.h ../include/libcmaes/esoptimizer.h ../include/libcmaes/esostrategy.h ../include/libcmaes/cmasolutions.h ../include/libcmaes/parameters.h ../include/libcmaes/cmaparameters.h ../include/libcmaes/cmastopcriteria.h ../include/libcmaes/ipopcmastrategy.h ../include/libcmaes/bipopcmastrategy.h ../include/libcmaes/covarianceupdate.h ../include/libcmaes/acovarianceupdate.h ../include/libcmaes/vdcmaupdate.h ../include/libcmaes/pwq_bound_strategy.h ../include/libcmaes/eigenmvn.h ../include/libcmaes/candidate.h ../include/libcmaes/genopheno.h ../include/libcmaes/noboundstrategy.h ../include/libcmaes/scaling.h ../include/libcmaes/llogging.h ../include/libcmaes/errstats.h ../include/libcmaes/pli.h ../include/libcmaes/contour.h

if HAVE_SURROG
libcmaes_la_SOURCES += surrcmaes.h surrogatestrategy.cc surrogatestrategy.h surrogates/rankingsvm.hpp surrogates/rsvm_surr_strategy.hpp
nobase_libcmaesinclude_HEADERS += ../include/libcmaes/surrcmaes.h ../include/libcmaes/surrogatestrategy.h ../include/libcmaes/surrogates/rankingsvm.hpp ../include/libcmaes/surrogates/rsvm_surr_strategy.hpp
endif

AM_CPPFLAGS=-I$(EIGEN3_INC) -I../include
AM_CXXFLAGS=-Wall -Wextra -g -O3
if !HAVE_CLANG
AM_CXXFLAGS += -fopenmp
endif

uninstall-hook:
	if test -d ${libcmaesincludedir}; then rm -rf ${libcmaesincludedir}; fi
