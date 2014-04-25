/**
 * CMA-ES, Covariance Matrix Evolution Strategy
 * Copyright (c) 2014 INRIA
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

#include "esoptimizer.h"
#include "cmastrategy.h"
#include <glog/logging.h>

using namespace libcmaes;

FitFunc cigtab = [](const double *x, const int N)
{
  int i;
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)
    sum += x[i]*x[i];
  return sum;
};

int main(int argc, char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr=1;
  google::SetLogDestination(google::INFO,"");
  //FLAGS_log_prefix=false;

  int dim = 5;
  int lambda = 10;
  CMAParameters<> cmaparams(dim,lambda);
  ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<>> cmaes(cigtab,cmaparams);
  cmaes.optimize();
}
