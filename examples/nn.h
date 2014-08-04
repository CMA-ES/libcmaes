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

#include "eo_matrix.h"
#include <vector>
#include <iostream>

double range = 1.0;
struct gen_rand {
  double factor;
public:
  gen_rand(double r=1.0) : factor(range/RAND_MAX) {}
  double operator()() {
    return rand() * factor;
  }
};

// build the network
//- input layer
//- hidden layer
//- top layer + loss function.

class nn_gradient
{
public:
  nn_gradient() {};
  ~nn_gradient() {};

  void initialize_gradients(const std::vector<int> &lsizes)
  {
    for (size_t i=0;i<lsizes.size()-1;i++)
      {
	int size_in = lsizes.at(i);
	int size_out = lsizes.at(i+1);
	dMat Gw = dMat::Zero(size_in,size_out);
	_GWs.push_back(Gw);
	dMat b = dMat::Zero(size_out,1);
	_Gbs.push_back(b);
      }
  }
  
  void accumulate(const std::vector<dMat> &GWupd,
		  const std::vector<dMat> &deltas)
  {
    int j = 0;
    for (int i=GWupd.size()-1;i>=0;i--)
      {
	_GWs.at(j) += GWupd.at(i).transpose();
	j++;
      }
    j = 0;
    for (int i=deltas.size()-1;i>=0;i--)
      {
	for (int k=0;k<deltas.at(i).cols();k++)
	  _Gbs.at(j) += deltas.at(i).col(k);
	j++;
      }
  }

  void grad_to_vec(const double &n,
		   const int &dim,
		   std::vector<double> &allgradient)
  {
    allgradient.reserve(dim);
    for (size_t i=0;i<_GWs.size();i++)
      {
	//std::cerr << "GW size=" << _GWs.at(i).size() << " / Gbs size=" << _Gbs.at(i).size() << std::endl;
	std::copy(_GWs.at(i).data(),_GWs.at(i).data()+_GWs.at(i).size(),std::back_inserter(allgradient));
	std::copy(_Gbs.at(i).data(),_Gbs.at(i).data()+_Gbs.at(i).size(),std::back_inserter(allgradient));
      }
    std::transform(allgradient.begin(),allgradient.end(),allgradient.begin(),
		   std::bind1st(std::multiplies<double>(),1.0/n));
    //std::cerr << "dim=" << dim << " / allgradient size=" << allgradient.size() << std::endl;
  }
  
  std::vector<dMat> _GWs;
  std::vector<dMat> _Gbs;
};

class nn
{
public:
  nn()
  {};
  
  nn(const std::vector<int> &lsizes,
     const int &unit=1, // 0: sigmoid, 1: tanh, 2: relu
     const bool &has_grad=false,
     const bool &r1=false,
     const bool &r2=false,
     const double &lambda1=0.1,
     const double &lambda2=0.1)
    :_lsizes(lsizes),_has_grad(has_grad),_unit(unit),_r1(r1),_r2(r2),_lambda1(lambda1),_lambda2(lambda2)
  {
    for (size_t i=0;i<_lsizes.size()-1;i++)
      {
	dMat weights = dMat::Random(_lsizes.at(i),_lsizes.at(i+1)) * sqrt(6.0/(_lsizes.at(i)+_lsizes.at(i+1)));
	if (_unit == 0)
	  weights *= 4.0;
	_lweights.push_back(weights);
	_lb.push_back(dVec::Zero(_lsizes.at(i+1)));
	_allparams_dim += _lweights.at(i).size() + _lsizes.at(i+1);
      }
    if (_has_grad)
      _grad.initialize_gradients(lsizes);
  }
  
  ~nn() {}

  dMat mtanh(const dMat &M)
  {
    dMat Mout = M;
    for (int i=0;i<Mout.rows();i++)
      for (int j=0;j<Mout.cols();j++)
	Mout(i,j) = tanh(M(i,j));
    return Mout;
  }

  dMat mtanh_derivative(const dMat &X)
  {
    dMat M = mtanh(X);
    M = M.cwiseProduct(M);
    return dMat::Constant(M.rows(),M.cols(),1)-M;
  }
  
  static dMat sigmoid(const dMat &M)
  {
    dMat expM(M.rows(),M.cols());
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    expM(i,j) = exp(-M(i,j));
	  }
      }
    dMat denom = expM + dMat::Constant(M.rows(),M.cols(),1);
    return denom.cwiseInverse();
  }

  dMat sigmoid_derivative(const dMat &X)
  {
    dMat M = sigmoid(X);
    return M.cwiseProduct((dMat::Constant(M.rows(),M.cols(),1)-M));
  }

  static dMat relu(const dMat &M)
  {
    dMat rM = M;
    return rM.cwiseMax(0);
  }

  dMat relu_derivative(const dMat &X)
  {
    dMat M = X;
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    M(i,j) = M(i,j) > 0 ? 1.0 : 0.0;
	  }
      }
    return M;
  }
  
  static dMat softmax(const dMat &M)
  {
    dMat maxM = M.colwise().maxCoeff();
    dMat expM(M.rows(),M.cols());
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    expM(i,j) = exp(M(i,j)-maxM(j));
	  }
      }
    dMat sums = expM.colwise().sum();
    // div by row vector.
    for (int i=0;i<expM.cols();i++)
      {
	for (int j=0;j<expM.rows();j++)
	  {
	    expM(j,i) /= sums(0,i);
	  }
      }
    return expM;
  }

  dMat get_loss(const dMat &prediction, const dMat &labels)
  {
    dMat M = labels.cwiseProduct(prediction);
    dMat Mc = M.colwise().sum();
    dMat logM(Mc.rows(),Mc.cols());
    for (int i=0;i<Mc.rows();i++)
      {
	for (int j=0;j<Mc.cols();j++)
	  {
	    logM(i,j) = Mc(i,j) == 0 ? -20.0 : log(Mc(i,j));
	  }
      }
    return (-1.0) * logM;
  }
    
  // forward pass over all features at once.
  void forward_pass(const dMat &features,
		    const dMat &labels)
  {
    to_matrices();
    for (size_t i=0;i<_lweights.size();i++)
      {
	dMat activation;
	if (i == 0)
	  activation = (_lweights.at(i).transpose() * features).colwise() + _lb.at(i);
	else activation = (_lweights.at(i).transpose() * _lfeatures).colwise() + _lb.at(i);
	if (_has_grad)
	  _activations.push_back(activation);
	if (i != _lweights.size()-1)
	  {
	    if (_unit == 1)
	      _lfeatures = mtanh(activation);
	    else if (_unit == 0)
	      _lfeatures = sigmoid(activation);
	    else if (_unit == 2)
	      _lfeatures = relu(activation);
	  }
	else _lfeatures = softmax(activation);
	if (_has_grad)
	  _predicts.push_back(_lfeatures);
      }
    
    // loss.
    if (labels.size() > 0) // training mode.
      {
	if (_has_grad)
	  {
	    dMat delta = _lfeatures - labels;
	    _deltas.push_back(delta);
	  }
	_loss = get_loss(_lfeatures,labels).mean();
	if (_r1)
	  {
	    double l1 = 0.0;
	    for (size_t i=0;i<_lweights.size();i++)
	      l1 += _lweights.at(i).cwiseAbs().sum();
	    _loss += _lambda1 * l1;
	  }
	if (_r2)
	  {
	    double l2 = 0.0;
	    for (size_t i=0;i<_lweights.size();i++)
	      l2 += _lweights.at(i).norm();
	    _loss += _lambda2 * l2;
	  }
      }
  }

  void back_propagate(const dMat &features)
  {
    std::vector<dMat> GWupd;
    int j = 0;
    for (int i=(int)_lweights.size()-1;i>=0;i--)
      {
	if (i > 0) // back propagate until hidden layer.
	  {
	    dMat delta = _lweights.at(i) * _deltas.back();
	    dMat der;
	    if (_unit == 0)
	      der = sigmoid_derivative(_activations.at(i-1));
	    else if (_unit == 1)
	      der = mtanh_derivative(_activations.at(i-1));
	    else if (_unit == 2)
	      der = relu_derivative(_activations.at(i-1));
	    delta = der.cwiseProduct(delta);
	    _deltas.push_back(delta);
	  }

	dMat GW;
	if (i == 0)
	  {
	    GW = _deltas.at(j) * features.transpose();
	  }
	else
	  {
	    GW = _deltas.at(j) * _predicts.at(i-1).transpose();
	  }  
	j++;
	GWupd.push_back(GW);
      }
    _grad.accumulate(GWupd,_deltas);
  }

  int sgd(const dMat &batchfeatures, const dMat &batchlabels,
	  const double &lrt)
  {
    // get loss.
    forward_pass(batchfeatures,batchlabels);
    
    // back propagate.
    back_propagate(batchfeatures);
    
    // update weights.
    grad_to_vec(batchfeatures.cols());
    std::vector<double> reg(_allparams.size());
    if (_r2)
      std::transform(_allparams.begin(),_allparams.end(),reg.begin(),std::bind1st(std::multiplies<double>(),-lrt));
    std::transform(_allgradient.begin(),_allgradient.end(),_allgradient.begin(),std::bind1st(std::multiplies<double>(),-lrt));
    std::transform(_allparams.begin(),_allparams.end(),_allgradient.begin(),_allparams.begin(),std::plus<double>());
    if (_r2)
      std::transform(_allparams.begin(),_allparams.end(),reg.begin(),_allparams.begin(),std::plus<double>());
    
    return 0;
  }
  
  void to_array()
  {
    _allparams.clear();
    _allparams.reserve(_allparams_dim);
    for (size_t i=0;i<_lweights.size();i++)
      {
	std::copy(_lweights.at(i).data(),_lweights.at(i).data()+_lweights.at(i).size(),std::back_inserter(_allparams));
	std::copy(_lb.at(i).data(),_lb.at(i).data()+_lsizes.at(i+1),std::back_inserter(_allparams));
      }
  }

  void to_matrices()
  {
    auto vit = _allparams.begin();
    for (size_t i=0;i<_lweights.size();i++)
      {
	std::copy(vit,vit+_lweights.at(i).size(),_lweights.at(i).data());
        vit += _lweights.at(i).size();
	std::copy(vit,vit+_lsizes.at(i+1),_lb.at(i).data());
	vit += _lsizes.at(i+1);
      }
  }

  dVec grad_to_vec(const double &n)
  {
    _allgradient.clear();
    _grad.grad_to_vec(n,_allparams.size(),_allgradient);
    dVec grad(_allparams.size());
    std::copy(_allgradient.begin(),_allgradient.end(),grad.data());
    return grad;
  }
  
  void clear_grad()
  {
    _grad = nn_gradient();
    _grad.initialize_gradients(_lsizes);
    _activations.clear();
    _deltas.clear();
    _predicts.clear();
  }

  bool grad_check(const dMat &gfeatures, const dMat &glabels)
  {
    double epsilon = 1e-4;
    int attempts = 1;
    while (attempts > 0)
      {
	_allparams.reserve(_allparams_dim);
	std::generate_n(std::back_inserter(_allparams), _allparams_dim, gen_rand());
	forward_pass(gfeatures,glabels);
	double cost = _loss;
	back_propagate(gfeatures);
	grad_to_vec(gfeatures.cols()); // in allgradient

	std::vector<double> allparams = _allparams;
	_has_grad = false;
	std::vector<double> numerical_gradient;
	for (int i=0;i<(int)_allparams_dim;i++)
	  {
	    std::vector<double> e(_allparams_dim,0.0);
	    e.at(i) = 2.0*epsilon;
	    std::vector<double> y1;
	    std::transform(allparams.begin(),allparams.end(),e.begin(),std::back_inserter(y1),std::plus<double>());
	    _allparams = y1;
	    //std::cerr << "old weight=" << allparams.at(i) << " / new weight=" << y1.at(i) << std::endl;
	    forward_pass(gfeatures,glabels);
	    double cost_epsilon = _loss;
	    double res = (cost_epsilon - cost) / (2.0*epsilon);
	    //std::cerr << "#i=" << i << " / cost=" << cost << " / cost_epsilon=" << cost_epsilon << " / res=" << res << std::endl;
	    numerical_gradient.push_back(res);
	  }

	double diff=0,total_diff=0,max_diff=-1.0;
	for (int i=0;i<(int)_allparams_dim;i++)
	  {
	    diff = fabs(numerical_gradient.at(i)-_allgradient.at(i));
	    std::cout << "i: " << i << " -- diff: " << diff
		      << " -- numerical gradient: " << numerical_gradient.at(i)
		      << " -- returned gradient: " << _allgradient.at(i) << std::endl;
	    max_diff = std::max(max_diff,diff);
	    total_diff += diff;
	  }
	std::cout << "attempt #" << attempts << " -- total diff: " << total_diff
		  << " -- max diff: " << max_diff << std::endl;
	
	if (max_diff > epsilon)
	  {
	    std::cout << "gradient cal failed: max diff too high\n";
	    return false;
	  }
	
	--attempts;
      }
    return true;
  }
  
  std::vector<int> _lsizes; /**< layer sizes. */
  std::vector<dMat> _lweights; /**< weight matrice, per layer. */
  std::vector<dVec> _lb;
  unsigned int _allparams_dim = 0;
  std::vector<double> _allparams; /**< all parameters, flat representation. */
  dMat _lfeatures;
  double _loss = std::numeric_limits<double>::max(); /**< current loss. */
  nn_gradient _grad;
  bool _has_grad = false;
  int _unit = false; /**< 0: sigmoid, 1: tanh, 2: relu. */
  std::vector<double> _allgradient; /**< flat representation. */
  std::vector<dMat> _activations; // only for bp.
  std::vector<dMat> _deltas;
  std::vector<dMat> _predicts;
  bool _r1 = false; /**< whether to use L1 regularization. */
  bool _r2 = false; /**< whether to use L2 regularization. */
  double _lambda1 = 0.1; /**< L1 factor. */
  double _lambda2 = 0.1; /**< L2 factor. */
};
