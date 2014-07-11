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
  double sigma = 0.1;
  std::vector<double> x0(dim,1.0);
  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
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
  double sigma = 0.1;
  std::vector<double> x0(dim,1.0);
  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  cmaparams._quiet = true;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  CMASolutions cmaksols = errstats<>::optimize_reduced_pk(fsphere,cmaparams,cmasols,6,0.1);
  std::cout << "iter: " << cmaksols._niter << std::endl;
  std::cout << "run status: " << cmaksols._run_status << std::endl;
  ASSERT_EQ(TOLHISTFUN,cmaksols._run_status);
  ASSERT_EQ(9,cmaksols.best_candidate()._x.size());
  std::cout << "fvalue: " << cmaksols.best_candidate()._fvalue << std::endl;
  std::cout << "x: " << cmaksols.best_candidate()._x.transpose() << std::endl;
}

TEST(optimize,optimize_fixed_p)
{
  FitFunc fsphere = [](const double *x, const int N)
    {
      double val = 0.0;
      for (int i=0;i<N;i++)
	val += x[i]*x[i];
      return val;
    };
  int dim = 10;
  double sigma = 0.1;
  std::vector<double> x0(dim,1.0);
  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  cmaparams._quiet = true;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  /*cmaparams.set_fixed_p(6,0.1);
    CMASolutions cmaksols = cmaes<>(fsphere,cmaparams);*/
  CMASolutions cmaksols = errstats<>::optimize_pk(fsphere,cmaparams,cmasols,6,0.1);
  std::cout << "iter: " << cmaksols._niter << std::endl;
  std::cout << "run status: " << cmaksols._run_status << std::endl;
  ASSERT_EQ(TOLHISTFUN,cmaksols._run_status);
  ASSERT_EQ(10,cmaksols.best_candidate()._x.size());
  std::cout << "fvalue: " << cmaksols.best_candidate()._fvalue << std::endl;
  std::cout << "x: " << cmaksols.best_candidate()._x.transpose() << std::endl;
}

TEST(pl,profile_likelihood_nocurve)
{
   FitFunc fsphere = [](const double *x, const int N)
    {
      double val = 0.0;
      for (int i=0;i<N;i++)
	val += x[i]*x[i];
      return val;
    };
   int dim = 10;
   double sigma = 0.1;
   std::vector<double> x0(dim,1.0);
   CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  cmaparams._quiet = true;
  cmaparams._seed = 1234;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  int k = 6;
  double fup = 0.1;
  int samplesize = 20;
  pli le = errstats<>::profile_likelihood(fsphere,cmaparams,cmasols,k,false,samplesize,fup);
  std::cout << "le fvalue: " << le._fvaluem.transpose() << std::endl;
  std::cout << "le xm: " << le._xm << std::endl;
  EXPECT_FLOAT_EQ(-0.31640971,le._min);
  EXPECT_FLOAT_EQ(0.31640971,le._max);
}

TEST(pl,profile_likelihood_curve)
{
   FitFunc fsphere = [](const double *x, const int N)
    {
      double val = 0.0;
      for (int i=0;i<N;i++)
	val += x[i]*x[i];
      return val;
    };
   int dim = 10;
   double sigma = 0.1;
  std::vector<double> x0(dim,1.0);
  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  cmaparams._quiet = true;
  cmaparams._seed = 4321;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  int k = 6;
  double fup = 0.1;
  int samplesize = 20;
  pli le = errstats<>::profile_likelihood(fsphere,cmaparams,cmasols,k,true,samplesize,fup);
  std::cout << "le fvalue: " << le._fvaluem.transpose() << std::endl;
  std::cout << "le xm: " << le._xm << std::endl;
  int mini, maxi;
  std::pair<double,double> mm = le.getMinMax(0.1,mini,maxi);
  EXPECT_FLOAT_EQ(-0.3313162,mm.first);
  EXPECT_FLOAT_EQ(0.3313162,mm.second);
}
