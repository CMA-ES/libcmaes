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
#include <iostream>

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

class customCMAStrategy : public CMAStrategy<CovarianceUpdate>
{
public:
  customCMAStrategy(FitFunc &func,
		    CMAParameters<> &parameters)
    :CMAStrategy<CovarianceUpdate>(func,parameters)
  {
  }

  ~customCMAStrategy() {}

  dMat ask()
  {
    return CMAStrategy<CovarianceUpdate>::ask();
  }

  void eval(const dMat &candidates,
	    const dMat &phenocandidates=dMat(0,0))
  {
    // custom eval.
    for (int r=0;r<candidates.cols();r++)
      {
	_solutions.get_candidate(r).set_x(candidates.col(r));
	if (phenocandidates.size()) // if candidates in phenotype space are given
	  _solutions.get_candidate(r).set_fvalue(_func(phenocandidates.col(r).data(),candidates.rows()));
	else _solutions.get_candidate(r).set_fvalue(_func(candidates.col(r).data(),candidates.rows()));
	
	//std::cerr << "candidate x: " << _solutions.get_candidate(r).get_x_dvec().transpose() << std::endl;
      }
    update_fevals(candidates.cols());
  }
  
  void tell()
  {
    CMAStrategy<CovarianceUpdate>::tell();
  }

  bool stop()
  {
    return CMAStrategy<CovarianceUpdate>::stop();
  }
  
};

int main(int argc, char *argv[])
{
  int dim = 10; // problem dimensions.
  std::vector<double> x0(dim,10.0);
  double sigma = 0.1;

  CMAParameters<> cmaparams(x0,sigma);
  //ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<>> optim(fsphere,cmaparams);
  ESOptimizer<customCMAStrategy,CMAParameters<>> optim(fsphere,cmaparams);
  
  while(!optim.stop())
    {
      dMat candidates = optim.ask();
      optim.eval(candidates);
      optim.tell();
      optim.inc_iter(); // important step: signals next iteration.
    }
  std::cout << optim.get_solutions() << std::endl;
}
