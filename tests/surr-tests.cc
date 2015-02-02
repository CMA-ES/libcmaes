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

#include "surrcmaes.h"
#include "surrogatestrategy.h"
#include <cmath>
#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

DEFINE_string(fname,"fsphere","name of the function to optimize");
DEFINE_string(dims,"10","comma-separated list of problem dimension");
//DEFINE_bool(tpa,false,"whether to use two-point adapation for step-size update");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_bool(no_quiet,false,"no intermediate output");
DEFINE_int32(runs,10,"number of runs for each configuration");
DEFINE_bool(no_exploit,false,"whether to exploit the surrogate model");
//DEFINE_int32(l,-1,"training set size (number of points)");
DEFINE_int32(lambdaprime,-1,"true objective function calls per iteration");
DEFINE_int32(prelambda,500,"number of pre-screened offprings sampled at every iteration");
DEFINE_int32(rsvm_iter,5e6,"number of iterations for optimizing the ranking SVM");

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
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

FitFunc elli = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=1;i<=N;i++)
    val += std::pow(1e6,(i-1)/(N-1))*x[i-1]*x[i-1];
  return val;
};

FitFunc schwefel = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    {
      double ival = 0.0;
      for (int j=0;j<i;j++)
	ival += x[j];
      val += ival*ival;
    }
  return val;
};

FitFunc schwefel14 = [](const double *x, const int N)
{
  return std::pow(schwefel(x,N),0.25);
};

FitFunc ackley = [](const double *x, const int N)
{
  double val = 0.0, isum = 0.0, icos = 0.0;
  for (int i=0;i<N;i++)
    {
      isum += x[i]*x[i];
      icos += std::cos(2.0*M_PI*x[i]);
    }
  val = -20.0*std::exp(-0.2*std::sqrt((1.0/static_cast<double>(N)) * isum)) - std::exp((1.0/static_cast<double>(N)) * icos) + 20.0 + std::exp(1);
  return val;
};

FitFunc rastrigin = [](const double *x, const int N)
{
  static double A = 10.0;
  double val = A*N;
  for (int i=0;i<N;i++)
    val += x[i]*x[i] - A*cos(2*M_PI*x[i]);
  return val;
};

std::map<std::string,FitFunc> mfuncs;
std::map<std::string,GradFunc> mgfuncs;

void tokenize(const std::string &str,
	      std::vector<std::string> &tokens,
	      const std::string &delim)
{
  
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delim, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delim, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delim, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delim, lastPos);
    }
}

template<template <class U, class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
  void set_optim_options(ESOptimizer<RSVMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> &optim,
			 const int &dim, const int &lambdaprime)
{
  if (FLAGS_no_exploit)
    optim.set_exploit(!FLAGS_no_exploit);
  optim.set_prelambda(FLAGS_prelambda);
  optim.set_lambdaprime(lambdaprime);
  //optim._rsvm_iter = FLAGS_rsvm_iter;
  optim._rsvm_iter = 50000*std::sqrt(dim);
  if (FLAGS_fname == "elli"
      || FLAGS_fname == "rosenbrock")
    optim.set_l(std::floor(70.0*std::sqrt(dim)));
    else if (FLAGS_fname == "rastrigin")
    {
      optim.set_theta_sel0(0.6);
      optim.set_theta_sel1(0.6);
    }
}

void run(const int &dim, const int &tpa, const std::string &alg,
	 double &fevals_avg, double &stddev, double &succ_runs)
{
  double x0min=-std::numeric_limits<double>::max(),x0max=-std::numeric_limits<double>::max();
  double sigma0 = -1;
  int lambda = -1;
  int lambdaprime = FLAGS_lambdaprime;
  if (dim == 2)
    {
      lambda = 6;
      lambdaprime = 3;
    }
  else if (dim == 4 || dim == 5)
    {
      lambda = 8;
      lambdaprime = 3;
    }
  else if (dim == 8 || dim == 10)
    {
      lambda = 10;
      lambdaprime = 3;
    }
  else if (dim == 16)
    {
      lambda = 11;
      lambdaprime = 3;
    }
  else if (dim == 20)
    {
      lambda = 12;
      lambdaprime = 4;
    }
  else if (dim == 32)
    {
      lambda = 14;
      lambdaprime = 4;
    }
  else if (dim == 40)
    {
      lambda = 15;
      lambdaprime = 5;
    }
  if (FLAGS_fname == "elli")
    {
      x0min = 1.0;
      x0max = 5.0;
      sigma0 = 2.0;
    }
  else if (FLAGS_fname == "rosenbrock")
    {
      x0min = -5.0;
      x0max = 5.0;
      sigma0 = 0.5;
    }
  else if (FLAGS_fname == "schwefel"
	   || FLAGS_fname == "schwefel14")
    {
      x0min = -10.0;
      x0max = 10.0;
      sigma0 = 10;
    }
  else if (FLAGS_fname == "ackley")
    {
      x0min = 1.0;
      x0max = 30.0;
      sigma0 = 14.5;
    }
  else if (FLAGS_fname == "rastrigin")
    {
      x0min = 5.0;
      x0max = 5.0;
      sigma0 = 2.0;
      if (dim == 2)
	{
	  lambda = 50;
	  lambdaprime = 25;
	}
      else if (dim == 5)
	{
	  lambda = 140;
	  lambdaprime = 70;
	}
    }
  fevals_avg = 0.0;
  succ_runs = 0.0;
  std::vector<double> vfevals;
  for (int r=0;r<FLAGS_runs;r++)
    {
      std::vector<double> x0(dim,-std::numeric_limits<double>::max());
      CMAParameters<> cmaparams(x0,sigma0,lambda);
      if (x0min != std::numeric_limits<double>::max())
	{
	  cmaparams.set_x0(x0min,x0max);
	}
      cmaparams.set_str_algo(alg);
      cmaparams.set_ftarget(1e-9);
      cmaparams.set_tpa(tpa);
      cmaparams.set_stopping_criteria(STAGNATION,false);
      cmaparams.set_stopping_criteria(TOLX,false);
      cmaparams.set_stopping_criteria(CONDITIONCOV,false);
      cmaparams.set_quiet(!FLAGS_no_quiet);
      ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(mfuncs[FLAGS_fname],cmaparams);
      set_optim_options(optim,dim,lambdaprime);
      optim.optimize();
      CMASolutions cmasols = optim.get_solutions();
      if (cmasols.best_candidate().get_fvalue() <= 1e-8)
	{
	  succ_runs++;
	  vfevals.push_back(cmasols.fevals());
	  fevals_avg += cmasols.fevals();
	}
    }
  fevals_avg /= succ_runs;
  for (double d: vfevals)
    stddev += std::pow(fevals_avg-d,2);
  stddev = std::sqrt((1.0/static_cast<double>(vfevals.size()))*stddev);
}

int main(int argc, char *argv[])
{
  mfuncs["fsphere"]=fsphere;
  mfuncs["rosenbrock"]=rosenbrock;
  mfuncs["elli"]=elli;
  mfuncs["schwefel"]=schwefel;
  mfuncs["schwefel14"]=schwefel14;
  mfuncs["ackley"]=ackley;
  mfuncs["rastrigin"]=rastrigin;

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // dims
  std::vector<std::string> vdims_str;
  tokenize(FLAGS_dims,vdims_str,",");
  std::vector<int> vdims;
  for (size_t i=0;i<vdims_str.size();i++)
    vdims.push_back(atoi(vdims_str.at(i).c_str()));

  // runs
  std::cout << "D\tfevals_avg\n";
  for (size_t d=0;d<vdims.size();d++)
    {
      int dim = vdims.at(d);
      double fevals_avg,stddev,succ_runs;
      run(dim,false,FLAGS_alg,fevals_avg,stddev,succ_runs);
      std::cout << dim << "\t" << fevals_avg <<  " +/- " << stddev << " (" << succ_runs << ")\n";
    }
}
