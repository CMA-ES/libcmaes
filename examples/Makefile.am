bin_PROGRAMS=sample_code sample_code_genopheno sample_code_pfunc sample_code_ask_tell sample_code_bounds sample_code_gradient sample_code_lscaling lorentzpeakbench sample_code_pffunc
sample_code_SOURCES=sample-code.cc
sample_code_genopheno_SOURCES=sample-code-genopheno.cc
sample_code_pfunc_SOURCES=sample-code-pfunc.cc
sample_code_ask_tell_SOURCES=sample-code-ask-tell.cc
lorentzpeakbench_SOURCES=lorentzpeakbench.cc
sample_code_bounds_SOURCES=sample-code-bounds.cc
sample_code_gradient_SOURCES=sample-code-gradient.cc
sample_code_lscaling_SOURCES=sample-code-lscaling.cc
sample_code_pffunc_SOURCES=sample-code-pffunc.cc

AM_CPPFLAGS=-I../src/ -I$(EIGEN3_INC)
AM_CXXFLAGS=-g -Wall -fopenmp
AM_LDFLAGS=-L../src
LDADD=-lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
