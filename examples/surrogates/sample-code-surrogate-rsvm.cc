/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "surrogatestrategy.h"
#include "opti_err.h"
#include "surrcmaes.h"
#include <map>
#include <iostream>

#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

FitFunc elli = [](const double *x, const int N)
{
  if (N == 1)
    return x[0] * x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += exp(log(1e3)*2.0*static_cast<double>(i)/static_cast<double>((N-1))) * x[i]*x[i];
  return val;
};

FitFunc rosenbrock = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N-1;i++)
    {
      val += 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
    }
  return val;
};

std::map<std::string,FitFunc> mfuncs;

DEFINE_string(fname,"fsphere","name of the function to optimize");
DEFINE_int32(dim,2,"problem dimension");
DEFINE_int32(lambda,-1,"number of offsprings");
DEFINE_int32(max_iter,-1,"maximum number of iteration (-1 for unlimited)");
DEFINE_int32(max_fevals,-1,"maximum budget as number of function evaluations (-1 for unlimited)");
DEFINE_double(sigma0,-1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_double(ftarget,-std::numeric_limits<double>::infinity(),"objective function target when known");
DEFINE_string(fplot,"","file where to store data for later plotting of results and internal states");
DEFINE_double(x0,-std::numeric_limits<double>::max(),"initial value for all components of the mean vector (-DBL_MAX for automated value)");
DEFINE_bool(noisy,false,"whether the objective function is noisy, automatically fits certain parameters");
DEFINE_bool(no_exploit,false,"whether to exploit the surrogate model");
DEFINE_int32(l,-1,"training set size (number of points)");
DEFINE_int32(lambdaprime,-1,"true objective function calls per iteration");
DEFINE_int32(prelambda,500,"number of pre-screened offprings sampled at every iteration");
DEFINE_int32(rsvm_iter,1e6,"number of iterations for optimizing the ranking SVM");

template<template <class U, class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
  void set_optim_options(ESOptimizer<RSVMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> &optim)
{
  if (FLAGS_no_exploit)
    optim.set_exploit(!FLAGS_no_exploit);
  optim.set_prelambda(FLAGS_prelambda);
  if (FLAGS_lambdaprime > 0)
    optim.set_lambdaprime(FLAGS_lambdaprime);
  optim._rsvm_iter = FLAGS_rsvm_iter;
  if (FLAGS_l > 0)
    optim.set_l(FLAGS_l);
}

int main(int argc, char *argv[])
{
  mfuncs["fsphere"]=fsphere;
  mfuncs["elli"]=elli;
  mfuncs["rosenbrock"]=rosenbrock;
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<double> x0(FLAGS_dim,FLAGS_x0);
  
  CMAParameters<> cmaparams(x0,FLAGS_sigma0);
  cmaparams.set_quiet(false);
  cmaparams.set_ftarget(FLAGS_ftarget);
  cmaparams.set_str_algo(FLAGS_alg);
  cmaparams.set_fplot(FLAGS_fplot);
  cmaparams.set_max_iter(FLAGS_max_iter);
  cmaparams.set_max_fevals(FLAGS_max_fevals);
  if (FLAGS_noisy)
    cmaparams.set_noisy();
  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(mfuncs[FLAGS_fname],cmaparams);
  set_optim_options(optim);
  optim.optimize();
  CMASolutions cmasols = optim.get_solutions();
  //CMASolutions cmasols = surrcmaes<>(mfuncs[FLAGS_fname],cmaparams); // can be used directly instead of defining ESOptimizer
  std::cout << cmasols << std::endl;
}
