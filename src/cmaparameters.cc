
#include "cmaparameters.h"
#include <math.h>
#include <iostream>
#include <glog/logging.h>

namespace libcmaes
{

  CMAParameters::CMAParameters(const int &dim, const int &lambda,
			       const int &max_iter)
    :Parameters(dim,lambda,max_iter)
  {
    _mu = ceil(_lambda / 2.0);
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

    //TODO: fast_ps, etc...
    
    _chi = sqrt(static_cast<double>(_dim))*(1.0-1.0/(4.0*_dim) + 1.0/(21.0*_dim*_dim));
  }

  CMAParameters::~CMAParameters()
  {
  }
  
}
