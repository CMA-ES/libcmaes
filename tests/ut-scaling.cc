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

#include "scaling.h"
#include "genopheno.h"
#include "esoptimizer.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace libcmaes;

TEST(linScalingStrategy,scale)
{
  int dim = 10;
  dVec scaling = dVec::Constant(dim,1.0);
  scaling(3) = 1e-3;
  dVec shift = dVec::Zero(dim);
  linScalingStrategy lsc(scaling,shift);
  dVec y = dVec::Constant(dim,1.0);
  dVec x;
  lsc.scale_to_internal(x,y);
  ASSERT_EQ(1e-3,x[3]);
  ASSERT_EQ(1.0,x[4]);
  lsc.scale_to_f(x,y);
  ASSERT_EQ(1.0,y[3]);
}

TEST(linScalingStrategy,compute)
{
  int dim = 10;
  double lbounds[dim],ubounds[dim];
  for (int i=0;i<dim;i++)
    {
      lbounds[i] = -1.0;
      ubounds[i] = 1.0;
    }
  linScalingStrategy lsc(lbounds,ubounds,dim);

  /*std::cerr << "scaling=" << lsc._scaling.transpose() << std::endl;
    std::cerr << "shift=" << lsc._shift.transpose() << std::endl;*/
  
  dVec y = dVec::Constant(dim,1.0);
  y[3] = -1.0;
  dVec x;
  lsc.scale_to_internal(x,y);
  ASSERT_EQ(10.0,x[4]);
  ASSERT_EQ(0.0,x[3]);
  lsc.scale_to_f(x,y);
  ASSERT_EQ(1.0,y[4]);
  ASSERT_EQ(-1.0,y[3]);

  y = dVec::Random(dim);
  x = dVec::Zero(dim);
  dVec yr;
  lsc.scale_to_internal(x,y);
  lsc.scale_to_f(x,yr);

  std::cerr << "y=" << y.transpose() << std::endl;
  std::cerr << "yr=" << yr.transpose() << std::endl;

  for (int i=0;i<dim;i++)
    ASSERT_FLOAT_EQ(y[i],yr[i]);
}

TEST(linScalingStrategy,bounds)
{
  int dim = 10;
  double lbounds[dim],ubounds[dim];
  for (int i=0;i<dim;i++)
    {
      lbounds[i] = -10.0;
      ubounds[i] = 10.0;
    }
  linScalingStrategy lsc(lbounds,ubounds,dim);

  for (int i=0;i<10;i++)
    {
      dVec y = dVec::Random(dim) * 10;
      dVec x,yr;
      lsc.scale_to_internal(x,y);
      lsc.scale_to_f(x,yr);
      
      /*std::cerr << "y=" << y.transpose() << std::endl;
	std::cerr << "yr=" << yr.transpose() << std::endl;*/
      
      for (int i=0;i<dim;i++)
	ASSERT_FLOAT_EQ(y[i],yr[i]);
    }
}

TEST(linScalingStrategy,bounds_overflow)
{
  int dim = 3;
  std::vector<double> lbounds = {-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),0.0};
  std::vector<double> ubounds = {std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),10.3818};
  linScalingStrategy lsc(&lbounds.front(),&ubounds.front(),dim);

  for (int i=0;i<10;i++)
    {
      dVec y = dVec::Random(dim) * 10;
      dVec x,yr;
      lsc.scale_to_internal(x,y);
      lsc.scale_to_f(x,yr);
      
      std::cerr << "y=" << y.transpose() << std::endl;
      std::cerr << "yr=" << yr.transpose() << std::endl;
	
      for (int i=0;i<dim;i++)
	ASSERT_FLOAT_EQ(y[i],yr[i]);
    }
}

TEST(linScalingStrategy,with_pwqbounds)
{
  // dummy genotype / phenotype transform functions.
  TransFunc genof = [](const double *ext, double *in, const int &dim)
    {
      for (int i=0;i<dim;i++)
	in[i] = 2.0*ext[i];
    };
  
  TransFunc phenof = [](const double *in, double *ext, const int &dim)
    {
      for (int i=0;i<dim;i++)
	ext[i] = 0.5*in[i];
    };

  for (int j=0;j<100;j++)
    {
      int dim = 3;
      std::vector<double> lbounds = {-2.0,-3.0,-4.0};
      std::vector<double> ubounds = {1.0,2.0,3.0};
      GenoPheno<pwqBoundStrategy,linScalingStrategy> gp(genof,phenof,&lbounds.front(),&ubounds.front(),dim);
      dVec candidates = dVec::Random(dim) + dVec::Constant(dim,1.0); // in [0,2]
      //std::cout << "candidates=" << candidates.transpose() << std::endl;
      dVec pcand = gp.pheno(candidates);
      dVec gcand = gp.geno(pcand);
      //std::cout << "pcand=" << pcand.transpose() << std::endl;
      //std::cout << "gcand=" << gcand.transpose() << std::endl;
      for (int i=0;i<dim;i++)
	ASSERT_FLOAT_EQ(gcand[i],candidates[i]);
    }
}

TEST(linScalingStrategy,vectorized_with_pwqbounds)
{
   // dummy genotype / phenotype transform functions.
  TransFunc genof = [](const double *ext, double *in, const int &dim)
    {
      for (int i=0;i<dim;i++)
	in[i] = 2.0*ext[i];
    };
  
  TransFunc phenof = [](const double *in, double *ext, const int &dim)
    {
      for (int i=0;i<dim;i++)
	ext[i] = 0.5*in[i];
    };

  for (int j=0;j<100;j++)
    {
      int dim = 3;
      int lambda = 10;
      std::vector<double> lbounds = {-2.0,-3.0,-4.0};
      std::vector<double> ubounds = {1.0,2.0,3.0};
      GenoPheno<pwqBoundStrategy,linScalingStrategy> gp(genof,phenof,&lbounds.front(),&ubounds.front(),dim);
      dMat candidates = dMat::Random(dim,lambda) + dMat::Constant(dim,lambda,1.0); // in [0,2]
      //std::cout << "candidates=" << candidates << std::endl;
      dMat pcand = gp.pheno(candidates);
      dMat gcand = gp.geno(pcand);
      //std::cout << "pcand=" << pcand << std::endl;
      //std::cout << "gcand=" << gcand << std::endl;
      for (int i=0;i<dim;i++)
	for (int k=0;k<lambda;k++)
	  ASSERT_FLOAT_EQ(gcand(i,k),candidates(i,k));
    }
}
