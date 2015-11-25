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
#include <cmath>
#include <fstream>
#include <iostream>
#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

double background(const double *x, const double *par)
{
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double lorentzianpeak(const double *x, const double *par)
{
  return (0.5*par[0]*par[1]/M_PI) /
    std::max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) + .25*par[1]*par[1]);
}

double blfunc(const double *x, const double *par)
{
  return background(x,par) + lorentzianpeak(x,&par[3]);
}


double points[201];
double values[201];

// chi2
FitFunc ff = [](const double *x, const int N)
{
  double sum = 0.0;
  for (int i=0;i<201;i++)
    {
      //std::cout << "x=" << points[i] << " / f=" << values[i] << std::endl;
      if (values[i] != 0.0)
	sum += pow((values[i]-blfunc(&points[i],x))/sqrt(values[i]),2);
    }
  //std::cout << "sum=" << sum << std::endl;
  return sum;
};

void loaddata(const std::string filename)
{
  std::ifstream fin(filename);
  std::string line;
  std::getline(fin,line); // skip header line.
  int a, b, i = 0;
  double c;
  while(fin >> a >> b >> c)
    {
      points[i] = c;
      values[i] = b;
      //std::cout << "x=" << points[i] << " / f=" << values[i] << std::endl;
      i++;
    }
  fin.close();
}

DEFINE_int32(lambda,-1,"number of offsprings");
DEFINE_double(sigma0,-1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_uint64(seed,0,"seed for random generator");
DEFINE_bool(le,false,"compute profile likelihood (confidence intervals)");
DEFINE_double(le_fup,0.1,"deviation from the minimum as the size of the confidence interval for profile likelihood computation");
DEFINE_double(le_delta,0.1,"tolerance factor around the fup confidence interval for profile likelihood computation");
DEFINE_int32(le_samplesize,10,"max number of steps of linesearch for computing the profile likelihood in every direction");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop");

int main(int argc, char *argv[])
{
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  double x0[6] = {1,1,1,6,.03,1};
  int dim = 6;
  loaddata("lorentzpeakbench.dat");
  CMAParameters<> cmaparams(dim,x0,FLAGS_sigma0,FLAGS_lambda,FLAGS_seed);
  /*if (FLAGS_x0 != std::numeric_limits<double>::min())
    cmaparams.set_x0(FLAGS_x0);*/
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
  else
    {
      std::cout << "unknown algorithm flavor " << FLAGS_alg << std::endl;
      exit(-1);
    }
  CMASolutions cmasols = cmaes<>(ff,cmaparams);
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";

  if (cmasols.run_status() >= 0 && FLAGS_le)
    {
      std::cerr << "Now computing confidence intervals around minimum for a deviation of " << FLAGS_le_fup << " (fval=" << cmasols.best_candidate().get_fvalue() + FLAGS_le_fup << ")\n";
      for (int k=0;k<dim;k++)
	errstats<>::profile_likelihood(ff,cmaparams,cmasols,k,false,
				       FLAGS_le_samplesize,FLAGS_le_fup,FLAGS_le_delta);
    }
  std::cout << "best solution: " << cmasols << std::endl;
  return cmasols.run_status();
}
