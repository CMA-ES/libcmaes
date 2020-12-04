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

#include <libcmaes/cmaes.h>
#include <iostream>

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

int main(int argc, char *argv[])
{
  const int dim = 10; // problem dimensions.
  std::vector<double> x0(dim,1.0);
  double sigma = 0.1;
  double lbounds[dim],ubounds[dim]; // arrays for lower and upper parameter bounds, respectively
  for (int i=0;i<dim;i++)
    {
      lbounds[i] = -2.0;
      ubounds[i] = 2.0;
    }
  GenoPheno<pwqBoundStrategy,linScalingStrategy> gp(lbounds,ubounds,dim);
  CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>> cmaparams(x0,sigma,-1,0,gp); // -1 for automatically decided lambda.
  cmaparams.set_algo(aCMAES);
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy,linScalingStrategy>>(fsphere,cmaparams);
  std::cout << "best solution: ";
  cmasols.print(std::cout,0,gp);
  std::cout << std::endl;
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  return cmasols.run_status();
}
