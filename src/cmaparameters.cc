
#include "cmaparameters.h"
#include <math.h>
#include <iostream>
#include <glog/logging.h>

namespace libcmaes
{

  CMAParameters::CMAParameters(const int &dim, const int &lambda,
			       const int &max_iter, const int &max_fevals,
			       const std::string &fplot,
			       const double &sigma_init, const double &x0,
			       const uint64_t &seed)
    :Parameters(dim,lambda,max_iter,max_fevals,x0,fplot,seed),_sigma_init(sigma_init),_nrestarts(9),_lazy_update(false),_lazy_value(0),_cm(1.0),_alphacov(2.0),_alphaminusold(0.5),_lambdamintarget(0.66)
  {
    _mu = floor(_lambda / 2.0);
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

    _csigma = (_muw+2.0)/(_dim+_muw+5.0);
    _cc = (4.0+_muw/static_cast<double>(_dim))/(_dim+4.0+2.0*_muw/static_cast<double>(_dim));
    _c1 = 2.0/(pow((_dim+1.3),2)+_muw);
    _cmu = std::min(1.0-_c1,2.0*(_muw-2.0+1.0/_muw)/(pow(_dim+2.0,2)+_muw));
    
    _dsigma = 1.0+_csigma+2.0*std::max(0.0,sqrt((_muw-1)/(_dim+1))-1);

    // constants used in covariance update.
    _fact_ps = sqrt(_csigma*(2.0-_csigma)*_muw);
    _fact_pc = sqrt(_cc * (2.0 - _cc) * _muw);
    
    _chi = sqrt(static_cast<double>(_dim))*(1.0-1.0/(4.0*_dim) + 1.0/(21.0*_dim*_dim));

    _lazy_value = 1.0/(_c1+_cmu)/_dim/10.0;

    // active cma.
    _deltamaxsigma = _alphaminusmin = std::numeric_limits<double>::max(); 
  }

  CMAParameters::~CMAParameters()
  {
    
  }
  
}
