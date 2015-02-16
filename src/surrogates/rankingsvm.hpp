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

/**
 * This is implementing a fast Ranking SVM algorithm, following:
 *        'Surrogate-Assisted Evolutionary Algorithms', Ilya Loshchilov, PhD Thesis, Universite Paris-Sud 11, 2013.
 *        http://www.loshchilov.com/phd.html
 *        see Chapter 4, p. 100.
 */

#ifndef RANKINGSVM_H
#define RANKINGSVM_H

#include "eo_matrix.h"
#include <vector>
#include <limits>
#include <cstdlib>
#include <random>
#include <iostream>

/**
 * \brief Kernel base class
 */
class SVMKernel
{
 public:
  SVMKernel() {};
  ~SVMKernel() {};

  double K(const dVec &x1, const dVec &x2);
  void init(const dMat &x) {};
};

/**
 * \brief linear kernel
 */
class LinearKernel : public SVMKernel
{
public:
  LinearKernel()
    :SVMKernel()
  {}

  ~LinearKernel() {}

  double K(const dVec &x1, const dVec &x2) { return x1.transpose()*x2; }
};

/**
 * \brief Polynomial kernel
 */
template <int d,int c=1>
class PolyKernel : public SVMKernel
{
public:
  PolyKernel()
  :SVMKernel()
  {}

  ~PolyKernel() {}

  double K(const dVec &x1, const dVec &x2) { return pow((x1.transpose()*x2 + c),d); }
};

/**
 * \brief Radial Basis Function kernel
 */
class RBFKernel : public SVMKernel
{
 public:
  RBFKernel()
    :SVMKernel()
    {}
  
  ~RBFKernel() {}

  double K(const dVec &x1, const dVec &x2) { return exp(-_gamma*((x1-x2).squaredNorm())); }

  void init(const dMat &x)
  {
    double avgdist = 0.0;
    for (int i=0;i<x.cols();i++)
      for (int j=i+1;j<x.cols();j++)
	avgdist += (x.col(i)-x.col(j)).norm();
    avgdist /= 0.5*(x.cols()*(x.cols()-1.0));
    double sigma = _sigma_a * std::pow(avgdist,_sigma_pow);
    _gamma = 1.0/(2.0*sigma*sigma);
    
    //debug
    //std::cout << "avgdist=" << avgdist << " / sigma=" << sigma << " / gamma=" << _gamma << std::endl;
    //debug    
  }
  
  double _gamma = 1.0;
  double _sigma_a = 1.0;
  double _sigma_pow = 1.0;
};

/**
 * \brief Ranking SVM algorithm with support for custom kernels
 */
template<class TKernel=RBFKernel>
class RankingSVM
{
 public:
  RankingSVM()
  {    
    _udist = std::uniform_real_distribution<>(0,1);
  }
  
  ~RankingSVM()
  {
  }

  /**
   * \brief trains a ranker from a set of points
   * @param x matrix in which every column represents a point from the training set
   * @param covinv the inverse sqrt covariance of the points in the training set, if
   *        available and needed, encode() function.
   * @param training set mean distribution if available along with covariance matrix,
   *        and needed, see encode() function
   * @see encode
   */
  void train(dMat &x,
	     const int &niter,
	     const dMat &covinv,
	     const dVec &xmean)
  {
    //debug
    //std::cout << "Learning RSVM with niter=" << niter << std::endl;
    //debug
    
    // init structures.
    int nalphas = x.cols()-1;
    _C = dMat::Constant(nalphas,1,_Cval);
    for (int i=0;i<nalphas;i++)
      _C(nalphas-1-i) = _Cval*pow(nalphas-i,2);
    _dKij = dMat::Zero(nalphas,nalphas);
    _alpha = dVec::Zero(nalphas);
    
    if (_encode)
      encode(x,covinv,xmean);
    compute_training_kernel(x);
    optimize(x,niter);
    
    //debug
    //std::cout << "alpha=" << _alpha.transpose() << std::endl;
    //debug
  }

  /**
   * \brief predicts a ranking from a learnt ranker
   * @param fit the final ranking fitted by the ranker
   * @param x_test points to be ranked, one per column of the matrix
   * @param x_train the initial training set, possibly used by kernel computation
   * @param covinv the inverse sqrt covariance of the points in the training set, if
   *        available and needed, encode() function.
   * @param training set mean distribution if available along with covariance matrix,
   *        and needed, see encode() function
   * @see encode
   */
  void predict(dVec &fit,
	       dMat &x_test,
	       dMat &x_train,
	       const dMat &covinv,
	       const dVec &xmean)
  {
    if (_alpha.size() == 0)
      return; // model is not yet trained.
    fit = dVec::Zero(x_test.cols());
    if (_encode)
      {
	encode(x_train,covinv,xmean);
	encode(x_test,covinv,xmean);
      }
#pragma omp parallel for
    for (int i=0;i<x_test.cols();i++)
      {
	dVec Kvals(x_train.cols());
	for (int j=0;j<x_train.cols();j++)
	  Kvals(j) = _kernel.K(x_test.col(i),x_train.col(j));
	double curfit = 0.0;
	for (int j=0;j<x_train.cols()-1;j++)
	  curfit += _alpha(j) * (Kvals(j)-Kvals(j+1));
	fit(i) = curfit;
      }
    
    //debug
    //std::cout << "fit=" << fit.transpose() << std::endl;
    //debug
  }

  /**
   * \brief encoding a set of point in a transformed space
   * @param x the points to be transformed, one per column of the matrix
   * @param covinv the inverse sqrt covariance of the points in the training set
   * @param training set mean distribution if available along with covariance matrix
   */
  void encode(dMat &x,
	      const dMat &covinv,
	      const dVec &xmean)
  {
    for (int i=0;i<x.cols();i++)
      x.col(i) -= xmean;
    if (covinv.cols() > 1)
      x = covinv * x;
    else
      {
	for (int i=0;i<x.cols();i++)
	  x.col(i) = covinv.cwiseProduct(x.col(i));
      }
  }

  /**
   * \brief pre-computation of the kernel values for every examples and coordinates
   * @param x training set as a point per column of the matrix
   */
  void compute_training_kernel(dMat &x)
  {
    _kernel.init(x);
    _K = dMat::Zero(x.cols(),x.cols());
#pragma omp parallel for
      for (int i=0;i<_K.rows();i++)
	for (int j=i;j<_K.cols();j++)
	  _K(i,j)=_K(j,i)=_kernel.K(x.col(i),x.col(j));
  
    //debug
    //std::cout << "K=" << _K << std::endl;
    //debug
  }

  /**
   * \brief optimizes a ranker's model given a training set x
   * @param training set as a set of points in column of the matrix
   * @param niter the number of iterations allowed for optimization
   */
  void optimize(const dMat &x,
		const int &niter)
  {
    // initialization of temporary variables
    dVec sum_alphas = dVec::Zero(_dKij.cols());
    dMat div_dKij = dMat::Zero(_dKij.rows(),_dKij.cols());
#pragma omp parallel
    {
#pragma omp for
      for (int i=0;i<_dKij.rows();i++)
	{
	  for (int j=0;j<_dKij.cols();j++)
	    {
	      _dKij(i,j) = _K(i,j) - _K(i,j+1) - _K(i+1,j) + _K(i+1,j+1);
	    }
	  double fact = _udist(_rng);
	  _alpha(i) = _C(i) * (0.95 + 0.05*fact);
	}
#pragma omp for
      for (int i=0;i<_dKij.rows();i++)
	{
	  double sum_alpha = 0.0;
	  for (int j=0;j<_dKij.cols();j++)
	    {
	      sum_alpha += _alpha(j) * _dKij(i,j);
	      div_dKij(i,j) = _dKij(i,j) / _dKij(j,j);
	    }
	  sum_alphas(i) = (_epsilon - sum_alpha) / _dKij(i,i);
	}
    }
    
    // optimize for niter
    double L=0.0;
    int i1 = 0;
    double old_alpha = 0.0, new_alpha = 0.0, delta_alpha = 0.0;
    for (int i=0;i<niter;i++)
      {
	i1 = i % _dKij.cols();
	old_alpha = _alpha(i1);
	new_alpha = old_alpha + sum_alphas(i1);
	new_alpha = std::max(std::min(new_alpha,_C(i1)),0.0);
	delta_alpha = new_alpha - old_alpha;
	double dL = delta_alpha * _dKij(i1,i1) * (sum_alphas(i1) - 0.5*delta_alpha + _epsilon);
	if (dL > 0)
	  {
	    sum_alphas -= delta_alpha * div_dKij.row(i1);
	    _alpha(i1) = new_alpha;
	  }
	L += dL;
      }
  }

  /**
   * \brief computes the ranker's error over a dataset
   * @param x_test testing dataset, with one point per column of the matrix
   * @param x_train the initial training dataset, used by encode() if needed
   * @param ref_fit the reference ranking fit against which to test the prediction
   * @param covinv the inverse sqrt covariance of the points in the training set, if
   *        available and needed, encode() function.
   * @param training set mean distribution if available along with covariance matrix,
   *        and needed, see encode() function
   * @see encode
   */
  double error(dMat &x_test,
	       dMat &x_train,
	       const dVec &ref_fit,
	       const dMat &covinv,
	       const dVec &xmean)
  {
    dVec fit;
    predict(fit,x_test,x_train,covinv,xmean);
    if (fit.size() == 0)
      return 1.0;
    double err = 0.0;
    double sum = 0.0;
    for (int i=0;i<ref_fit.size();i++)
      {
	for (int j=0;j<ref_fit.size();j++)
	  {
	    if (i != j)
	      {
		err += ((ref_fit(i) > ref_fit(j) && fit(i) < fit(j)) || (ref_fit(i) < ref_fit(j) && fit(i) > fit(j))) ? 1 : 0;
		sum++;
	      }
	  }
      }
    err /= static_cast<double>((ref_fit.size()*ref_fit.size())-ref_fit.size());
    return err;
  }

 public:
  bool _encode = false; /**< whether to use encoding from inverse sqrt covariance matrix of points. */
  
  dMat _K; /**< pre-computed matrix of kernel values for a given training set. */
  dVec _alpha; /**< vector of Ranking SVM parameters over ranking constraints. */
  dMat _dKij;
  dMat _C; /**< constraint violation weights. */
  double _Cval = 1e6; /**< constraing violation base weight value. */
  double _epsilon = 1.0;
  
  TKernel _kernel; /**< kernel class. */

  std::mt19937 _rng;
  std::uniform_real_distribution<> _udist;
};

#endif
