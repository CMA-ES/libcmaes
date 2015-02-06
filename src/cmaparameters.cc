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

#include "cmaparameters.h"
#include <cmath>
#include <iostream>

namespace libcmaes
{
  template<class TGenoPheno>
  CMAParameters<TGenoPheno>::CMAParameters(const int &dim,
					   const double *x0,
					   const double &sigma,
					   const int &lambda,
					   const uint64_t &seed,
					   const TGenoPheno &gp)
    :Parameters<TGenoPheno>(dim,x0,lambda,seed,gp),_sigma_init(sigma),_nrestarts(9),_lazy_update(false),_lazy_value(0),_cm(1.0),_alphacov(2.0),_alphaminusold(0.5),_lambdamintarget(0.66),_alphaminusmin(1.0)
  {
    initialize_parameters();
  }

  template<class TGenoPheno>
  CMAParameters<TGenoPheno>::CMAParameters(const std::vector<double> &x0,
					   const double &sigma,
					   const int &lambda,
					   const uint64_t &seed,
					   const TGenoPheno &gp)
    :Parameters<TGenoPheno>(x0.size(),&x0.front(),lambda,seed,gp),_sigma_init(sigma),_nrestarts(9),_lazy_update(false),_lazy_value(0),_cm(1.0),_alphacov(2.0),_alphaminusold(0.5),_lambdamintarget(0.66),_alphaminusmin(1.0)
  {
    initialize_parameters();
  }

  template<class TGenoPheno>
  CMAParameters<TGenoPheno>::CMAParameters(const std::vector<double> &x0,
					   const std::vector<double> &sigma,
					   const int &lambda,
					   const std::vector<double> &lbounds,
					   const std::vector<double> &ubounds,
					   const uint64_t &seed)
    :Parameters<TGenoPheno>(x0.size(),&x0.front(),lambda,seed,TGenoPheno()),_nrestarts(9),_lazy_update(false),_lazy_value(0),_cm(1.0),_alphacov(2.0),_alphaminusold(0.5),_lambdamintarget(0.66),_alphaminusmin(1.0)
  {
    dVec scaling = dVec::Constant(x0.size(),1.0).cwiseQuotient(Eigen::Map<dVec>(const_cast<double*>(&sigma.front()),sigma.size()));
    dVec shift = dVec::Constant(x0.size(),0.0);
    TGenoPheno gp(scaling,shift,&lbounds.front(),&ubounds.front()); // XXX: is only effective when GenoPheno has linScalingStrategy
    this->set_gp(gp);
    _sigma_init = *std::min_element(sigma.begin(),sigma.end());
    initialize_parameters();
  }
  
  template <class TGenoPheno>
  CMAParameters<TGenoPheno>::~CMAParameters()
  {
  }

  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::initialize_parameters()
  {
    _mu = floor(Parameters<TGenoPheno>::_lambda / 2.0);
    _weights = dVec::Zero(_mu);
    double sum_weights = 0.0, sq_weights = 0.0;
    for (int i=0;i<_mu;i++)
      {
	_weights[i] = log(_mu+1)-log(i+1);
	sum_weights += _weights[i];
	sq_weights += _weights[i] * _weights[i];
      }
    _weights /= sum_weights;

    //debug
    //std::cerr << "weights: " << _weights << std::endl;
    //debug
    
    _muw = sum_weights*sum_weights / sq_weights;

    if (Parameters<TGenoPheno>::_dim < 1000)
      _csigma = (_muw+2.0)/(Parameters<TGenoPheno>::_dim+_muw+5.0);
    else _csigma = (std::sqrt(_muw) + 2.0) / (std::sqrt(Parameters<TGenoPheno>::_dim) + std::sqrt(_muw) + 3.0);
    _cc = (4.0+_muw/static_cast<double>(Parameters<TGenoPheno>::_dim))/(Parameters<TGenoPheno>::_dim+4.0+2.0*_muw/static_cast<double>(Parameters<TGenoPheno>::_dim));
    
    _c1 = 2.0/(pow((Parameters<TGenoPheno>::_dim+1.3),2)+_muw);
    _cmu = 2.0*(_muw-2.0+1.0/_muw)/(pow(Parameters<TGenoPheno>::_dim+2.0,2)+_muw);
    _cmu = std::min(1.0-_c1,_cmu);
    
    _dsigma = 1.0+_csigma+2.0*std::max(0.0,sqrt((_muw-1)/(Parameters<TGenoPheno>::_dim+1))-1);
    
    // constants used in covariance update.
    _fact_ps = sqrt(_csigma*(2.0-_csigma)*_muw);
    _fact_pc = sqrt(_cc * (2.0 - _cc) * _muw);
    
    _chi = sqrt(static_cast<double>(Parameters<TGenoPheno>::_dim))*(1.0-1.0/(4.0*Parameters<TGenoPheno>::_dim) + 1.0/(21.0*Parameters<TGenoPheno>::_dim*Parameters<TGenoPheno>::_dim));
    
    _lazy_value = 1.0/(_c1+_cmu)/Parameters<TGenoPheno>::_dim/10.0;

    // active cma.
    _deltamaxsigma = std::numeric_limits<double>::max();

    // uncertainty handling.
    this->_rlambda = std::max(0.1,2.0/Parameters<TGenoPheno>::_lambda);
    this->_alphathuh = 1 + 2.0/(Parameters<TGenoPheno>::_dim+10.0);
  }
  
  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::set_noisy()
  {
    static double factor = 0.2;
    static double lfactor = 5.0; // lambda factor.
    Parameters<TGenoPheno>::_lambda *= lfactor;
    initialize_parameters(); // reinit parameters.
    _c1 *= factor;
    _cmu = std::min(1.0-_c1,2.0*factor*(_muw-2.0+1.0/_muw)/(pow(Parameters<TGenoPheno>::_dim+2.0,2)+_muw));
  }

  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::set_tpa(const int &b)
  {
    this->_tpa = b;
    if (this->_tpa == 2)
      _dsigma = 4 - 3.6 / std::sqrt(Parameters<TGenoPheno>::_dim);
    else if (!this->_vd)
      _dsigma = 1.0+_csigma+2.0*std::max(0.0,sqrt((_muw-1)/(Parameters<TGenoPheno>::_dim+1))-1); // default
    else _dsigma = 1.0+_csigma+2.0*std::max(0.0,sqrt((_muw-1)/(Parameters<TGenoPheno>::_dim+1))-1); // vd
  }
  
  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::reset_as_fixed(const int &k)
  {
    Parameters<TGenoPheno>::_dim--;
    removeElement(Parameters<TGenoPheno>::_x0min,k); // XXX: could go into parameters.cc
    removeElement(Parameters<TGenoPheno>::_x0max,k);
    removeElement(_weights,k);
  }

  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::set_sep()
  {
    _sep = true;
    _c1 *= (Parameters<TGenoPheno>::_dim+2.0)/3.0;
    _cmu = std::min(1.0-_c1,(0.3 + _muw - 2.0 + 1.0/_muw) / (this->_dim + 4.0*std::sqrt(this->_dim) + _muw/2.0));
    _lazy_value = 1.0/(_c1+_cmu)/Parameters<TGenoPheno>::_dim/10.0;
    _cc = (1.0+1.0/static_cast<double>(Parameters<TGenoPheno>::_dim) + _muw/static_cast<double>(Parameters<TGenoPheno>::_dim))/(std::sqrt(Parameters<TGenoPheno>::_dim)+1.0/static_cast<double>(Parameters<TGenoPheno>::_dim)+2.0*_muw/static_cast<double>(Parameters<TGenoPheno>::_dim));
  }
  
  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::set_vd()
  {
    if (this->_algo != 12 && this->_algo != 13 && this->_algo != 14)
      {
	std::cerr << "[Warning]: set_vd on non VD algorithm " << this->_algo << ". Not activating VD update\n";
	return;
      }
    _vd = true;
    _csigma = std::sqrt(_muw)/(2.0*(std::sqrt(Parameters<TGenoPheno>::_dim) + std::sqrt(_muw)));
    _c1 *= (Parameters<TGenoPheno>::_dim-5)/6.0;
    _cmu = std::min(1.0-_c1,(Parameters<TGenoPheno>::_dim-5)/6.0*(2.0*(_muw-2.0+1.0/_muw)/(pow(Parameters<TGenoPheno>::_dim+2.0,2)+_muw)));
    if (this->_tpa == 0)
      _dsigma = 1.0+_csigma+2.0*std::max(0.0,sqrt((_muw-1)/(Parameters<TGenoPheno>::_dim+1))-1);
    _fact_ps = sqrt(_csigma*(2.0-_csigma)*_muw);
  }

  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::set_fixed_p(const int &index, const double &value)
  {
    Parameters<TGenoPheno>::set_fixed_p(index,value);
    double ndim = Parameters<TGenoPheno>::_dim-Parameters<TGenoPheno>::_fixed_p.size();
    _chi = sqrt(ndim)*(1.0-1.0/(4.0*ndim) + 1.0/(21.0*ndim*ndim));
    _lazy_value = 1.0/(_c1+_cmu)/ndim/10.0;
  }

  template <class TGenoPheno>
  void CMAParameters<TGenoPheno>::unset_fixed_p(const int &index)
  {
    Parameters<TGenoPheno>::unset_fixed_p(index);
    double ndim = Parameters<TGenoPheno>::_dim-Parameters<TGenoPheno>::_fixed_p.size();
    _chi = sqrt(ndim)*(1.0-1.0/(4.0*ndim) + 1.0/(21.0*ndim*ndim));
    _lazy_value = 1.0/(_c1+_cmu)/ndim/10.0;
  }
  
  template class CMAParameters<GenoPheno<NoBoundStrategy>>;
  template class CMAParameters<GenoPheno<pwqBoundStrategy>>;
  template class CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
