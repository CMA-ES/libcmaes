bin_PROGRAMS=sample_code sample_code_genopheno sample_code_pfunc sample_code_ask_tell sample_code_ask_tell_uh sample_code_bounds sample_code_gradient sample_code_lscaling sample_code_lscaling_sigmas sample_code_pffunc
sample_code_SOURCES=sample-code.cc
sample_code_genopheno_SOURCES=sample-code-genopheno.cc
sample_code_pfunc_SOURCES=sample-code-pfunc.cc
sample_code_ask_tell_SOURCES=sample-code-ask-tell.cc
sample_code_ask_tell_uh_SOURCES=sample-code-ask-tell-uh.cc
sample_code_bounds_SOURCES=sample-code-bounds.cc
sample_code_gradient_SOURCES=sample-code-gradient.cc
sample_code_lscaling_SOURCES=sample-code-lscaling.cc
sample_code_lscaling_sigmas_SOURCES=sample-code-lscaling-sigmas.cc
sample_code_pffunc_SOURCES=sample-code-pffunc.cc

if HAVE_SURROG
bin_PROGRAMS += sample_code_surrogate1 test_rsvm
sample_code_surrogate1_SOURCES=surrogates/sample-code-surrogate1.cc
test_rsvm_SOURCES=surrogates/test-rsvm.cc
endif

if HAVE_GFLAGS
bin_PROGRAMS += lorentzpeakbench
lorentzpeakbench_SOURCES=lorentzpeakbench.cc
if HAVE_SURROG
bin_PROGRAMS += sample_code_simple_surrogate_rsvm sample_code_surrogate_rsvm
sample_code_simple_surrogate_rsvm_SOURCES=surrogates/sample-code-simple-surrogate-rsvm.cc
sample_code_surrogate_rsvm_SOURCES=surrogates/sample-code-surrogate-rsvm.cc
endif
endif

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC) $(GFLAGS_CFLAGS)
AM_CXXFLAGS=-g -Wall
if !HAVE_CLANG
AM_CXXFLAGS+=-fopenmp
endif
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
