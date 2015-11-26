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

#include "cmaes.h"
#include "errstats.h"
#include <map>
#include <random>
#include <limits>
#include <iostream>
#include <cmath>
#include "llogging.h"

//#define STRIP_FLAG_HELP 1
#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

#include <assert.h>

using namespace libcmaes;

std::random_device rd;
std::normal_distribution<double> norm(0.0,1.0);
std::cauchy_distribution<double> cauch(0.0,1.0);
std::mt19937 gen;
std::string boundtype = "none";

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

bool compEp(const double &a, const double &b, const double &epsilon)
{
  return fabs(a-b) <= epsilon;
}

dMat orthogonalBasis(const int N)
{
  static dMat b = dMat::Zero(N,N);
  static bool initialized = false;

  if (initialized)
    return b;
  
  std::random_device rd;
  std::normal_distribution<double> norm(0.0,1.0);
  std::mt19937 gen(rd()); //TODO: seed ?
  initialized = true;
  double sp = 0.0;
  
  for (int i=0;i<N;i++)
    {
      /* sample from Gaussian. */
      for (int j=0;j<N;j++)
	b(i,j) = norm(gen);
      /* substract projection of previous vectors */
      for (int j=i-1;j>=0;--j)
	{
	  sp = 0.0;
	  for (int k=0;k<N;++k)
	    sp += b(i,k)*b(j,k); // scalar product.
	  for (int k=0;k<N;k++)
	    b(i,k) -= sp * b(j,k);
	}
      sp = 0.0;
      for (int k=0;k<N;++k)
	sp += b(i,k)*b(i,k);
      for (int k=0;k<N;++k)
	b(i,k) /= sqrt(sp);
    }
  return b;
};

// frand for debug.
typedef std::mt19937 myrng; // Mersenne Twister.
myrng rng;
std::uniform_real_distribution<> uint_dist(0,1);
FitFunc frand = [](const double *x, const int N)
{
  return uint_dist(rng);
};

// classical test functions for single-objective optimization problems.
FitFunc ackleys = [](const double *x, const int N)
{
  return -20.0*exp(-0.2*sqrt(0.5*(x[0]*x[0]+x[1]*x[1]))) - exp(0.5*(cos(2.0*M_PI*x[0]) + cos(2.0*M_PI*x[1]))) + 20.0 + exp(1.0);
};

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

GradFunc grad_fsphere = [](const double *x, const int N)
{
  dVec grad(N);
  for (int i=0;i<N;i++)
    grad(i) = 2.0*x[i];
  return grad;
};

FitFunc cigtab = [](const double *x, const int N)
{
  int i;
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)
    sum += x[i]*x[i];
  return sum;
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

GradFunc grad_rosenbrock = [](const double *x, const int N)
{
  dVec grad = dVec::Zero(N);
  for (int i=0;i<N-1;i++)
    {
      grad(i) = -400.0*x[i]*(x[i+1]-x[i]*x[i])-2.0*(1.0-x[i]);
      grad(i+1) += 200.0*(x[i+1]-x[i]*x[i]);
    }
  return grad;
};

FitFunc beale = [](const double *x, const int N)
{
  return pow(1.5-x[0]+x[0]*x[1],2) + pow(2.25 - x[0] + x[0]*x[1]*x[1],2) + pow(2.625-x[0]+x[0]*pow(x[1],3),2);
};

FitFunc goldstein_price = [](const double *x, const int N)
{
  return (1.0 + pow(x[0] + x[1] + 1.0,2)*(19.0-14.0*x[0]+3*x[0]*x[0]-14.0*x[1]+6*x[0]*x[1]+3*x[1]*x[1]))*(30.0+pow(2.0*x[0]-3.0*x[1],2)*(18.0-32.0*x[0]+12.0*x[0]*x[0]+48.0*x[1]-36.0*x[0]*x[1]+27.0*x[1]*x[1]));
};

FitFunc booth = [](const double *x, const int N)
{
  return pow(x[0]+2.0*x[1]-7,2) + pow(2*x[0]+x[1]-5,2);
};

FitFunc bukin = [](const double *x, const int N)
{
  return 100.0 * sqrt(fabs(x[1]-0.01*x[0]*x[0])) + 0.01*fabs(x[0]+10.0);
};

FitFunc matyas = [](const double *x, const int N)
{
  return 0.26*(x[0]*x[0]+x[1]*x[1])-0.48*x[0]*x[1];
};

FitFunc levi = [](const double *x, const int N)
{
  return pow(sin(3*M_PI*x[0]),2) + pow(x[0]-1,2) * (1.0+pow(sin(3*M_PI*x[1]),2)) + pow(x[1]-1,2)*(1.0+pow(sin(2.0*M_PI*x[1]),2));
};

FitFunc camel = [](const double *x, const int N)
{
  return 2.0*x[0]*x[0] - 1.05*pow(x[0],4) + pow(x[0],6)/6.0 + x[0]*x[1] + x[1]*x[1];
};

FitFunc easom = [](const double *x, const int N)
{
  return -cos(x[0])*cos(x[1])*exp(-(pow((x[0]-M_PI),2)+pow((x[1]-M_PI),2)));
};

FitFunc crossintray = [](const double *x, const int N)
{
  return -0.0001*pow(fabs(sin(x[0])*sin(x[1])*exp(fabs(100.0-sqrt(x[0]*x[0]+x[1]*x[1])/M_PI)))+1.0,0.1);
};

FitFunc eggholder = [](const double *x, const int N)
{
  return -(x[1]+47)*sin(sqrt(fabs(x[1]+0.5*x[0]+47.0))) - x[0]*sin(sqrt(fabs(x[0]-(x[1] + 47.0))));
};

FitFunc holdertable = [](const double *x, const int N)
{
  return -fabs(sin(x[0])*cos(x[1])*exp(fabs(1.0-sqrt(x[0]*x[0]+x[1]*x[1])/M_PI)));
};

FitFunc mccormick = [](const double *x, const int N)
{
  return sin(x[0]+x[1])+pow(x[0]-x[1],2) - 1.5*x[0] + 2.5*x[1] + 1.0;
};

FitFunc schaffer2 = [](const double *x, const int N)
{
  return 0.5 + (pow(sin(x[0]*x[0]-x[1]*x[1]),2)-0.5) / pow(1.0+0.001*(x[0]*x[0]+x[1]*x[1]),2);
};

FitFunc schaffer4 = [](const double *x, const int N)
{
  return 0.5 + (cos(sin(fabs(x[0]*x[0]-x[1]*x[1])))-0.5) / pow(1.0+0.001*(x[0]*x[0]+x[1]*x[1]),2);
};

FitFunc styblinski_tang = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += pow(x[i],4) - 16.0*x[i]*x[i] + 5.0*x[i];
  return 0.5*val;
};

FitFunc rastrigin = [](const double *x, const int N)
{
  static double A = 10.0;
  double val = A*N;
  for (int i=0;i<N;i++)
    val += x[i]*x[i] - A*cos(2*M_PI*x[i]);
  return val;
};
std::vector<double> rastx0(10,-std::numeric_limits<double>::max()); // auto x0 in [-4,4].
CMAParameters<> rastrigin_params(rastx0,5.0,400,1234); // 1234 is seed.

FitFunc elli = [](const double *x, const int N)
{
  if (N == 1)
    return x[0] * x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += exp(log(1e3)*2.0*static_cast<double>(i)/static_cast<double>((N-1))) * x[i]*x[i];
  return val;
};

GradFunc grad_elli = [](const double *x, const int N)
{
  dVec grad(N);
  if (N == 1)
    {
      grad(0) = 2.0*x[0];
      return grad;
    }
  for (int i=0;i<N;i++)
    grad(i) = exp(log(1e3)*2.0*static_cast<double>(i)/static_cast<double>((N-1)))*2.0*x[i];
  return grad;
};

FitFunc tablet = [](const double *x, const int N)
{
  double val = 1e6*x[0]*x[0];
  for (int i=1;i<N;i++)
    val += x[i]*x[i];
  return val;
};

GradFunc grad_tablet = [](const double *x, const int N)
{
  dVec grad(N);
  grad(0) = 1e6*2.0*x[0];
  for (int i=0;i<N;i++)
    grad(i) = 2.0*x[i];
  return grad;
};

FitFunc cigar = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=1;i<N;i++)
    val += x[i]*x[i];
  val *= 1e6;
  val += x[0]*x[0];
  return val;
};

FitFunc ellirot = [](const double *x, const int N)
{
  if (N == 1)
    return x[0]*x[0];
  dMat b = orthogonalBasis(N);
  double val = 0.0;
  for (int i=0;i<N;i++)
    {
      double y = 0.0;
      for (int k=0;k<N;k++)
	y += b(i,k)*x[k];
      val += exp(log(1e3)*2.0*static_cast<double>(i)/(N-1)) * y*y;
    }
  return val;
};

FitFunc diffpow = [](const double *x, const int N)
{
  if (N == 1)
    return x[0]*x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += pow(fabs(x[i]),2.0+10.0*static_cast<double>(i)/(N-1));
  return val;
};

FitFunc diffpowrot = [](const double *x, const int N)
{
  if (N == 1)
    return x[0]*x[0];
  dMat b = orthogonalBasis(N);
  double val = 0.0;
  for (int i=0;i<N;i++)
    {
      double y = 0.0;
      for (int k=0;k<N;k++)
	y += b(i,k)*x[k];
      val += pow(fabs(y),2.0+10.0*static_cast<double>(i)/(N-1));
    }
  return val;
};

FitFunc hardcos = [](const double *x, const int N)
{
  double sum = 0.0;
  for (int i=0;i<N;i++)
    sum += x[i]*x[i];
  sum*=(cos(sum)+2.0);
  return sum;
};

// uncertainty handling testing functions
FitFunc fsphere_uhc = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  val += cauch(gen); // noise
  return val;
};

FitFunc elli_uh = [](const double *x, const int N)
{
  if (N == 1)
    return x[0] * x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += pow(10,6*static_cast<double>(i)/static_cast<double>((N-1))) * x[i]*x[i];
  val += norm(gen); // noise
  return val;
};

FitFunc elli_uhc = [](const double *x, const int N)
{
  if (N == 1)
    return x[0] * x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += pow(10,6*static_cast<double>(i)/static_cast<double>((N-1))) * x[i]*x[i];
  val += cauch(gen); // noise
  return val;
};

std::map<std::string,FitFunc> mfuncs;
std::map<std::string,GradFunc> mgfuncs;
std::map<std::string,Candidate> msols;
std::map<std::string,CMAParameters<>> mparams;
std::map<std::string,FitFunc>::const_iterator mit;
std::map<std::string,Candidate>::const_iterator fmit;
std::map<std::string,CMAParameters<>>::const_iterator pmit;

void fillupfuncs()
{
  mfuncs["frand"]=frand;
  mfuncs["ackleys"]=ackleys;
  msols["ackleys"]=Candidate(0.0,dVec::Constant(2,0));
  mfuncs["fsphere"]=fsphere;
  mgfuncs["fsphere"]=grad_fsphere;
  msols["fsphere"]=Candidate(0.0,dVec::Constant(20,0));
  mfuncs["cigtab"]=cigtab;
  mfuncs["rosenbrock"]=rosenbrock;
  mgfuncs["rosenbrock"]=grad_rosenbrock;
  msols["rosenbrock"]=Candidate(0.0,dVec::Constant(20,1));
  mfuncs["beale"]=beale;
  mfuncs["goldstein_price"]=goldstein_price;
  mfuncs["booth"]=booth;
  mfuncs["bukin"]=bukin;
  mfuncs["matyas"]=matyas;
  mfuncs["levi"]=levi;
  mfuncs["camel"]=camel;
  mfuncs["easom"]=easom;
  mfuncs["crossintray"]=crossintray;
  mfuncs["eggholder"]=eggholder;
  mfuncs["holdertable"]=holdertable;
  mfuncs["mccormick"]=mccormick;
  mfuncs["schaffer2"]=schaffer2;
  mfuncs["schaffer4"]=schaffer4;
  mfuncs["styblinski_tang"]=styblinski_tang;
  mfuncs["rastrigin"]=rastrigin;
  msols["rastrigin"]=Candidate(0.0,dVec::Constant(10,1));
  rastrigin_params.set_x0(5.0);
  mparams["rastrigin"]=rastrigin_params;
  mfuncs["elli"]=elli;
  mgfuncs["elli"]=grad_elli;
  msols["elli"]=Candidate(0.0,dVec::Constant(10,0));
  mfuncs["tablet"]=tablet;
  mgfuncs["tablet"]=grad_tablet;
  msols["tablet"]=Candidate(0.0,dVec::Constant(10,0));
  mfuncs["cigar"]=cigar;
  msols["cigar"]=Candidate(0.0,dVec::Constant(10,0));
  mfuncs["ellirot"]=ellirot;
  msols["ellirot"]=Candidate(0.0,dVec::Constant(10,0));
  mfuncs["diffpow"]=diffpow;
  mfuncs["diffpowrot"]=diffpowrot;
  mfuncs["hardcos"]=hardcos;
  mfuncs["fsphere_uhc"]=fsphere_uhc;
  mfuncs["elli_uh"]=elli_uh;
  mfuncs["elli_uhc"]=elli_uhc;
}

void printAvailFuncs()
{
  std::cout << "available functions: ";
  for (auto imap: mfuncs)
    std::cout << imap.first << " ";
  std::cout << std::endl;
}

// command line options.
DEFINE_string(fname,"fsphere","name of the function to optimize");
DEFINE_int32(dim,2,"problem dimension");
DEFINE_int32(lambda,-1,"number of offsprings");
DEFINE_int32(max_iter,-1,"maximum number of iteration (-1 for unlimited)");
DEFINE_int32(max_fevals,-1,"maximum budget as number of function evaluations (-1 for unlimited)");
DEFINE_bool(list,false,"returns a list of available functions");
DEFINE_bool(all,false,"test on all functions");
DEFINE_double(epsilon,1e-10,"epsilon on function result testing, with --all");
DEFINE_string(fplot,"","file where to store data for later plotting of results and internal states");
DEFINE_bool(full_fplot,false,"whether to activate full legacy plot");
DEFINE_double(sigma0,-1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_double(x0,-std::numeric_limits<double>::max(),"initial value for all components of the mean vector (-DBL_MAX for automated value)");
DEFINE_uint64(seed,0,"seed for random generator");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_bool(lazy_update,false,"covariance lazy update");
//DEFINE_string(boundtype,"none","treatment applied to bounds, none or pwq (piecewise linear / quadratic) transformation");
DEFINE_double(lbound,std::numeric_limits<double>::max()/-1e2,"lower bound to parameter vector");
DEFINE_double(ubound,std::numeric_limits<double>::max()/1e2,"upper bound to parameter vector");
DEFINE_bool(quiet,false,"no intermediate output");
DEFINE_bool(le,false,"whether to return profile likelihood error bounds around the minimum");
DEFINE_double(le_fup,0.1,"deviation from the minimum as the size of the confidence interval for profile likelihood computation");
DEFINE_double(le_delta,0.1,"tolerance factor around the fup confidence interval for profile likelihood computation");
DEFINE_int32(le_samplesize,10,"max number of steps of linesearch for computing the profile likelihood in every direction");
DEFINE_int32(le_maxiters,1e4,"max number of iterations in search for profile likelihood points");
DEFINE_bool(noisy,false,"whether the objective function is noisy, automatically fits certain parameters");
DEFINE_string(contour,"","two comma-separated variable indexes to which passes a contour to be computed as a set of additional points");
DEFINE_int32(contour_p,4,"number of contour points, must be >= 4");
DEFINE_double(contour_fup,0.1,"value from the minimum at which to find contour points");
DEFINE_bool(linscaling,false,"whether to automatically scale parameter space linearly so that parameter sensitivity is similar across all dimensions (requires -lbound and/or -ubound");
DEFINE_double(ftarget,-std::numeric_limits<double>::infinity(),"objective function target when known");
DEFINE_int32(restarts,9,"maximum number of restarts, applies to IPOP and BIPOP algorithms");
DEFINE_bool(with_gradient,false,"whether to use the function gradient when available in closed form");
DEFINE_bool(with_num_gradient,false,"whether to use numerical gradient injection");
DEFINE_bool(with_edm,false,"whether to compute expected distance to minimum when optimization has completed");
DEFINE_bool(with_stds,false,"whether to compute and print the standard deviation for every parameter");
DEFINE_bool(with_errors,false,"whether to compute the errors");
DEFINE_bool(with_corr,false,"whether to compute and print the correlation matrix (may not fit in memory in large-scale settings)");
DEFINE_bool(mt,false,"whether to use parallel evaluation of objective function");
DEFINE_bool(initial_fvalue,false,"whether to compute initial objective function value at x0");
DEFINE_int32(elitist,0,"whether to activate elistism, 0: deactivated, 1: reinjects best seen candidate, 2: initial elitism, reinjects x0, 3: on restart scheme, useful when optimizer appears to converge to a value that is higher than the best value reported along the way");
DEFINE_int32(max_hist,-1,"maximum stored history, helps mitigate the memory usage though preventing the 'stagnation' criteria to trigger");
DEFINE_bool(no_stagnation,false,"deactivate stagnation stopping criteria");
DEFINE_bool(no_tolx,false,"deactivate tolX stopping criteria");
DEFINE_bool(no_automaxiter,false,"deactivate automaxiter stopping criteria");
DEFINE_bool(no_tolupsigma,false,"deactivate tolupsigma stopping criteria");
DEFINE_bool(uh,false,"activate uncertainty handling of objective function");
DEFINE_int32(tpa,1,"whether to use two-point adapation for step-size update, 0: no, 1: auto, 2: yes");
DEFINE_double(tpa_dsigma,-1,"set two-point adaptation dsigma (use with care)");

template <class TGenoPheno=GenoPheno<NoBoundStrategy,NoScalingStrategy>>
CMASolutions cmaes_opt()
{
  std::vector<double> lbounds = {FLAGS_lbound},ubounds = {FLAGS_ubound};
  if (FLAGS_lbound != std::numeric_limits<double>::max()/-1e2 || FLAGS_ubound != std::numeric_limits<double>::max()/1e2)
    {
      boundtype = "pwq";
      lbounds = std::vector<double>(FLAGS_dim);
      ubounds = std::vector<double>(FLAGS_dim);
      for (int i=0;i<FLAGS_dim;i++)
	{
	  lbounds[i] = FLAGS_lbound;
	  ubounds[i] = FLAGS_ubound;
	}
    }
  TGenoPheno gp(&lbounds.at(0),&ubounds.at(0),FLAGS_dim);
  std::vector<double> x0(FLAGS_dim,FLAGS_x0);
  CMAParameters<TGenoPheno> cmaparams(x0,FLAGS_sigma0,FLAGS_lambda,FLAGS_seed,gp);
  cmaparams.set_max_iter(FLAGS_max_iter);
  cmaparams.set_max_fevals(FLAGS_max_fevals);
  cmaparams.set_restarts(FLAGS_restarts);
  cmaparams.set_fplot(FLAGS_fplot);
  cmaparams.set_full_fplot(FLAGS_full_fplot);
  cmaparams.set_lazy_update(FLAGS_lazy_update);
  cmaparams.set_quiet(FLAGS_quiet);
  cmaparams.set_tpa(FLAGS_tpa);
  cmaparams.set_gradient(FLAGS_with_gradient || FLAGS_with_num_gradient);
  cmaparams.set_edm(FLAGS_with_edm);
  cmaparams.set_mt_feval(FLAGS_mt);
  cmaparams.set_initial_fvalue(FLAGS_initial_fvalue);
  cmaparams.set_elitism(FLAGS_elitist);
  cmaparams.set_max_hist(FLAGS_max_hist);
  cmaparams.set_uh(FLAGS_uh);
  if (FLAGS_tpa_dsigma > 0.0)
    cmaparams.set_tpa_dsigma(FLAGS_tpa_dsigma);
  if (FLAGS_ftarget != -std::numeric_limits<double>::infinity())
    cmaparams.set_ftarget(FLAGS_ftarget);
  if (FLAGS_noisy)
    cmaparams.set_noisy();
  //cmaparams.set_str_algo(FLAGS_alg);
  if (FLAGS_alg == "cmaes")
    cmaparams.set_algo(CMAES_DEFAULT);
  else if (FLAGS_alg == "ipop")
    cmaparams.set_algo(IPOP_CMAES);
  else if (FLAGS_alg == "bipop")
    cmaparams.set_algo(BIPOP_CMAES);
  else if (FLAGS_alg == "acmaes")
    cmaparams.set_algo(aCMAES);
  else if (FLAGS_alg == "aipop")
    cmaparams.set_algo(aIPOP_CMAES);
  else if (FLAGS_alg == "abipop")
    cmaparams.set_algo(aBIPOP_CMAES);
  else if (FLAGS_alg == "sepcmaes")
    cmaparams.set_algo(sepCMAES);
  else if (FLAGS_alg == "sepipop")
    cmaparams.set_algo(sepIPOP_CMAES);
  else if (FLAGS_alg == "sepbipop")
    cmaparams.set_algo(sepBIPOP_CMAES);
  else if (FLAGS_alg == "sepacmaes")
    cmaparams.set_algo(sepaCMAES);
  else if (FLAGS_alg == "sepaipop")
    cmaparams.set_algo(sepaIPOP_CMAES);
  else if (FLAGS_alg == "sepabipop")
    cmaparams.set_algo(sepaBIPOP_CMAES);
  else if (FLAGS_alg == "vdcma")
    cmaparams.set_algo(VD_CMAES);
  else if (FLAGS_alg == "vdipopcma")
    cmaparams.set_algo(VD_IPOP_CMAES);
  else if (FLAGS_alg == "vdbipopcma")
    cmaparams.set_algo(VD_BIPOP_CMAES);
  else
    {
      LOG(ERROR) << "unknown algorithm flavor " << FLAGS_alg << std::endl;
      exit(-1);
    }
  if (FLAGS_no_automaxiter)
    cmaparams.set_stopping_criteria(AUTOMAXITER,false);
  if (FLAGS_no_stagnation)
    cmaparams.set_stopping_criteria(STAGNATION,false);
  if (FLAGS_no_tolx)
    cmaparams.set_stopping_criteria(TOLX,false);
  if (FLAGS_no_tolupsigma)
    cmaparams.set_stopping_criteria(TOLUPSIGMA,false);
  
  CMASolutions cmasols;
  GradFunc gfunc = nullptr;
  if (FLAGS_with_gradient)
    gfunc = mgfuncs[FLAGS_fname];
  PlotFunc<CMAParameters<TGenoPheno>,CMASolutions> pffunc = CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultFPFunc;
  if (cmaparams.dim() > 300)
    {
      pffunc = [](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
	{
	  std::string sep = " ";
	  fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << sqrt(cmasols.max_eigenv()/cmasols.min_eigenv()) << sep << cmasols.elapsed_last_iter() << std::endl;
	  return 0;
	};
    }
  cmasols = cmaes<>(mfuncs[FLAGS_fname],cmaparams,CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultPFunc,gfunc,cmasols,pffunc);
  std::cout << "Minimization completed in " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  if (cmasols.run_status() >= 0 && FLAGS_le)
    {
      std::cout << "Now computing confidence interval around minimum for a deviation of " << FLAGS_le_fup << " (fval=" << cmasols.best_candidate().get_fvalue() + FLAGS_le_fup << ")\n";
      for (int k=0;k<FLAGS_dim;k++)
	errstats<TGenoPheno>::profile_likelihood(mfuncs[FLAGS_fname],cmaparams,cmasols,k,false,
						 FLAGS_le_samplesize,FLAGS_le_fup,FLAGS_le_delta,FLAGS_le_maxiters);
    }
  if (!FLAGS_contour.empty())
    {
      cmaparams.set_quiet(true);
      std::vector<std::string> contour_indexes_str = split(FLAGS_contour,',');
      std::pair<int,int> contour_indexes;
      contour_indexes.first = atoi(contour_indexes_str.at(0).c_str());
      contour_indexes.second = atoi(contour_indexes_str.at(1).c_str());
      std::cout << "Now computing contour passing through point (" << contour_indexes.first << "," << contour_indexes.second << ")\n";
      contour ct = errstats<TGenoPheno>::contour_points(mfuncs[FLAGS_fname],contour_indexes.first,contour_indexes.second,
							FLAGS_contour_p,FLAGS_contour_fup,cmaparams,cmasols,FLAGS_le_delta,FLAGS_le_maxiters);
      std::cout << ct << std::endl;
    }
  std::cout << "Done!\n";
  if (cmasols.run_status() < 0)
    LOG(INFO) << "optimization failed with termination criteria " << cmasols.run_status() << " -- " << cmasols.status_msg() << std::endl;
  else LOG(INFO) << "optimization succeeded with termination criteria " << cmasols.run_status() << " -- " << cmasols.status_msg() << std::endl;
  LOG(INFO) << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  if (cmasols.best_candidate().get_x_size() <= 1000)
    {
      cmasols.print(std::cout,false,cmaparams.get_gp());
      std::cout << std::endl;
    }
  if (FLAGS_with_edm)
    LOG(INFO) << "EDM=" << cmasols.edm() << " / EDM/fm=" << cmasols.edm() / cmasols.best_candidate().get_fvalue() << std::endl;
  if (FLAGS_with_stds)
    std::cout << "stds=" << cmasols.stds(cmaparams).transpose() << std::endl;
  if (FLAGS_with_errors)
    std::cout << "errors=" << cmasols.errors(cmaparams).transpose() << std::endl;
  if (FLAGS_with_corr)
    std::cout << "correlation=" << cmasols.corr() << std::endl;
  return cmasols;
}

int main(int argc, char *argv[])
{
  gflags::ParseCommandLineFlags(&argc, &argv, true);
#ifdef HAVE_GLOG
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr=1;
  google::SetLogDestination(google::INFO,"");
  //FLAGS_log_prefix=false;
#endif
  
  fillupfuncs();
  
  gen = std::mt19937(rd());
  gen.seed(static_cast<uint64_t>(time(nullptr)));

  if (FLAGS_list)
    {
      printAvailFuncs();
      exit(1);
    }
  else if (FLAGS_all)
    {
      mit = mfuncs.begin();
      while(mit!=mfuncs.end())
	{
	  if ((fmit=msols.find((*mit).first))==msols.end())
	    {
	      ++mit;
	      continue;
	    }
	  int dim = msols[(*mit).first].get_x_dvec().rows();
	  std::vector<double> x0(dim,FLAGS_x0);
	  CMAParameters<> cmaparams(x0,FLAGS_sigma0,FLAGS_lambda);
	  cmaparams.set_max_iter(FLAGS_max_iter);
	  if ((pmit=mparams.find((*mit).first))!=mparams.end())
	    cmaparams = (*pmit).second;
	  cmaparams.set_quiet(true);
	  cmaparams.set_lazy_update(FLAGS_lazy_update);
	  if (FLAGS_alg == "cmaes")
	    cmaparams.set_algo(CMAES_DEFAULT);
	  else if (FLAGS_alg == "ipop")
	    cmaparams.set_algo(IPOP_CMAES);
	  else if (FLAGS_alg == "bipop")
	    cmaparams.set_algo(BIPOP_CMAES);
	  else if (FLAGS_alg == "acmaes")
	    cmaparams.set_algo(aCMAES);
	  else if (FLAGS_alg == "aipop")
	    cmaparams.set_algo(aIPOP_CMAES);
	  else if (FLAGS_alg == "abipop")
	    cmaparams.set_algo(aBIPOP_CMAES);
	  CMASolutions cmasols = cmaes<>(mfuncs[(*mit).first],cmaparams);
	  Candidate c = cmasols.best_candidate();
	  //TODO: check on solution in x space.
	  if (compEp(c.get_fvalue(),(*fmit).second.get_fvalue(),FLAGS_epsilon))
	    LOG(INFO) << (*mit).first << " -- OK\n";
	  else LOG(INFO) << (*mit).first << " -- FAILED - f-value=" << c.get_fvalue() << " / expected f-value=" << (*fmit).second.get_fvalue() << std::endl;
	  ++mit;
	}
      exit(1);
    }
  
  if ((mit=mfuncs.find(FLAGS_fname))==mfuncs.end())
    {
      LOG(ERROR) << FLAGS_fname << " function does not exist, run with --list to get the list of all functions. Exiting.\n";
      printAvailFuncs();
      exit(1);
    }
  CMASolutions cmasols;
  if (boundtype == "none")
    {
      if (!FLAGS_linscaling)
	cmasols = cmaes_opt<>();
      else cmasols = cmaes_opt<GenoPheno<NoBoundStrategy,linScalingStrategy>>();
    }
  else if (boundtype == "pwq")
    {
      if (!FLAGS_linscaling)
	cmasols = cmaes_opt<GenoPheno<pwqBoundStrategy>>();
      else cmasols = cmaes_opt<GenoPheno<pwqBoundStrategy,linScalingStrategy>>();
    }
  else
    {
      LOG(ERROR) << "Unknown boundtype " << boundtype << std::endl;
      exit(-1);
    }
}
