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
#include "surrogatestrategy.h"
#include <iostream>

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

CSurrFunc ftrain = [](const std::vector<Candidate> &c, const dMat &cov)
{
  //do nothing
  return 0;
};

SurrFunc fpredict = [](std::vector<Candidate> &c, const dMat &cov)
{
  // fill up with real fvalue.
  for (size_t i=0;i<c.size();i++)
    c.at(i).set_fvalue(fsphere(c.at(i).get_x_ptr(),c.at(i).get_x_size()));
  return 0;
};

int main(int argc, char *argv[])
{
  int dim = 10; // problem dimensions.
  std::vector<double> x0(dim,10.0);
  double sigma = 0.1;

  CMAParameters<> cmaparams(x0,sigma);
  ESOptimizer<SimpleSurrogateStrategy<CMAStrategy>,CMAParameters<>> optim(fsphere,cmaparams);
  optim.set_ftrain(ftrain);
  optim.set_fpredict(fpredict);
  optim.set_exploit(false); // test mode.
  
  while(!optim.stop())
    {
      dMat candidates = optim.ask();
      optim.eval(candidates);
      optim.tell();
      optim.inc_iter(); // important step: signals next iteration.
    }
  std::cout << optim.get_solutions() << std::endl;
}
