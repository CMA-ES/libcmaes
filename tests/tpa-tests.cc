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
#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

DEFINE_string(fname,"fsphere","name of the function to optimize");
DEFINE_string(dims,"10","comma-separated list of problem dimension");
//DEFINE_bool(tpa,false,"whether to use two-point adapation for step-size update");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_int32(runs,10,"number of runs for each configuration");
//DEFINE_bool(with_gradient,false,"whether to use the function gradient when available in closed form");

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

void run(const int &dim, const bool &gi, const int &tpa, const std::string &alg,
	 double &fevals_avg, double &succ_runs)
{
  fevals_avg = 0.0;
  succ_runs = 0.0;
  for (int r=0;r<FLAGS_runs;r++)
    {
      std::vector<double> x0(dim,-std::numeric_limits<double>::max());
      CMAParameters<> cmaparams(x0,-1);
      cmaparams.set_gradient(gi);
      cmaparams.set_str_algo(alg);
      cmaparams.set_ftarget(1e-8);
      cmaparams.set_tpa(tpa);
      cmaparams.set_stopping_criteria(STAGNATION,false);
      cmaparams.set_stopping_criteria(TOLX,false);
      //cmaparams.set_quiet(false);
      GradFunc gfunc = nullptr;
      if (gi)
	gfunc = mgfuncs[FLAGS_fname];
      CMASolutions cmasols = cmaes<>(mfuncs[FLAGS_fname],cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,gfunc);
      if (cmasols.best_candidate().get_fvalue() <= 1e-8)
	succ_runs++;
      fevals_avg += cmasols.fevals();
    }
  fevals_avg /= succ_runs;
}

int main(int argc, char *argv[])
{
  mfuncs["fsphere"]=fsphere;
  mgfuncs["fsphere"]=grad_fsphere;
  mfuncs["rosenbrock"]=rosenbrock;
  mgfuncs["rosenbrock"]=grad_rosenbrock;
  mfuncs["elli"]=elli;
  mgfuncs["elli"]=grad_elli;

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // dims
  std::vector<std::string> vdims_str;
  tokenize(FLAGS_dims,vdims_str,",");
  std::vector<int> vdims;
  for (size_t i=0;i<vdims_str.size();i++)
    vdims.push_back(atoi(vdims_str.at(i).c_str()));

  // dsigma
  /*std::vector<std::string> vds_str;
  tokenize(FLAGS_tpa_dsigmas,vds_str,",");
  std::vector<double> vds;
  for (size_t i=0;i<vds_str.size();i++)
  vds.push_back(strtod(vds_str.at(i).c_str(),NULL));*/

  // runs
  std::cout << "D\tfevals_avg\t\tfevals_avg_tpa\t\tfevals_avg_gi\t\tfevals_avg_gi_tpa\n";
  for (size_t d=0;d<vdims.size();d++)
    {
      int dim = vdims.at(d);
      double fevals_avg,fevals_avg_tpa,fevals_avg_gi,fevals_avg_gi_tpa;
      double succ_runs,succ_runs_tpa,succ_runs_gi,succ_runs_gi_tpa;
      run(dim,false,0,FLAGS_alg,fevals_avg,succ_runs);
      run(dim,false,2,FLAGS_alg,fevals_avg_tpa,succ_runs_tpa);
      run(dim,true,0,FLAGS_alg,fevals_avg_gi,succ_runs_gi);
      run(dim,true,2,FLAGS_alg,fevals_avg_gi_tpa,succ_runs_gi_tpa);
      std::cout << dim << "\t" << fevals_avg <<  " (" << succ_runs << ")\t\t" 
		<< fevals_avg_tpa << " (" << succ_runs_tpa << ")\t\t" 
		<< fevals_avg_gi << " (" << succ_runs_gi << ")\t\t" 
		<< fevals_avg_gi_tpa << "(" << succ_runs_gi_tpa << ")\n";
    }
}
