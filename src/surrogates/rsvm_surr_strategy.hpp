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
#include "opti_err.h"
#include "rankingsvm.hpp"

#ifndef RSVMSURROGATESTRATEGY_H
#define RSVMSURROGATESTRATEGY_H

namespace libcmaes
{

  // Interfacing candidates with an Eigen matrix of double x and a vector of objective function value.
  void to_mat_vec(std::vector<Candidate> &cp,
		  dMat &x, dVec &fvalues,
		  const bool &train)
  {
    if (train)
      std::sort(cp.begin(),cp.end(),
		[](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() > c2.get_fvalue();}); // descending sort
    x = dMat(cp.at(0).get_x_size(),cp.size());
    fvalues = dVec(cp.size());
    for (int i=0;i<(int)cp.size();i++)
      {
	x.col(i) = cp.at(i).get_x_dvec().transpose();
	fvalues(i) = cp.at(i).get_fvalue();
      }
  }

  template <class TGenoPheno> using eostrat = ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >;
  
  template<template <class U, class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class RSVMSurrogateStrategy : public ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>
    {
    public:
    RSVMSurrogateStrategy(FitFunc &func,
			  CMAParameters<TGenoPheno> &parameters)
          :ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>(func,parameters)
    {
      this->_train = [this](const std::vector<Candidate> &c, const dMat &cov)
      {
	if (c.empty())
	  return 0;
	dMat x;
	dVec fvalues;
	std::vector<Candidate> cp = c;
	to_mat_vec(cp,x,fvalues,true);
	dVec xmean = eostrat<TGenoPheno>::get_solutions().xmean();
	_rsvm = RankingSVM<RBFKernel>();
	_rsvm._encode = true;
	_rsvm.train(x,_rsvm_iter,cov,xmean);
	return 0;
      };
      this->_predict = [this](std::vector<Candidate> &c, const dMat &cov)
      {
	dMat x_test(c.at(0).get_x_size(),c.size());
	for (int i=0;i<(int)c.size();i++)
	  x_test.col(i) = c.at(i).get_x_dvec().transpose();
	
	dMat x_train;
	dVec fvalues;
	std::vector<Candidate> tset = this->_tset;
	to_mat_vec(tset,x_train,fvalues,true);
	
	dVec fit;
	dVec xmean = eostrat<TGenoPheno>::get_solutions().xmean();
	_rsvm.predict(fit,x_test,x_train,cov,xmean);
	if (fit.size() != 0)
	  for (int i=0;i<(int)c.size();i++)
	    c.at(i).set_fvalue(fit(i));
	return 0;
      };
    }
      
      ~RSVMSurrogateStrategy() {}
      
      RankingSVM<RBFKernel> _rsvm;
      int _rsvm_iter = 1e6; /**< number of iterations for optimizing the ranking SVM */
  };

}

#endif
