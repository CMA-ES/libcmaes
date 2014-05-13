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

#include "errstats.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace libcmaes;

TEST(rearrangecmasol,reset_as_fixed)
{
  FitFunc fsphere = [](const double *x, const int N)
    {
      double val = 0.0;
      for (int i=0;i<N;i++)
	val += x[i]*x[i];
      return val;
    };
  int dim = 10;
  CMAParameters<> cmaparams(dim);
  cmaparams._quiet = true;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  CMASolutions copsols = cmasols;
  copsols.reset_as_fixed(6);
  ASSERT_EQ(9,copsols._cov.rows());
  ASSERT_EQ(9,copsols._cov.cols());
  ASSERT_EQ(9,copsols._xmean.size());
  /*std::cout << cmasols._xmean.transpose() << std::endl;
    std::cout << copsols._xmean.transpose() << std::endl;*/
  CMAParameters<> copparams = cmaparams;
  cmaparams.reset_as_fixed(6);
  
}

TEST(optimize,optimize_pk)
{
  FitFunc fsphere = [](const double *x, const int N)
    {
      double val = 0.0;
      for (int i=0;i<N;i++)
	val += x[i]*x[i];
      return val;
    };
  int dim = 10;
  CMAParameters<> cmaparams(dim);
  cmaparams._quiet = false;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  CMASolutions cmaksols = errstats<>::optimize_pk(fsphere,cmaparams,cmasols,6,0.1);
  std::cout << "iter: " << cmaksols._niter << std::endl;
  std::cout << "run status: " << cmaksols._run_status << std::endl;
  ASSERT_EQ(TOLHISTFUN,cmaksols._run_status);
  ASSERT_EQ(9,cmaksols.best_candidate()._x.size());
  std::cout << "fvalue: " << cmaksols.best_candidate()._fvalue << std::endl;
  std::cout << "x: " << cmaksols.best_candidate()._x.transpose() << std::endl;
}
