
#include "cmasolutions.h"
#include <iostream>

namespace libcmaes
{

  CMASolutions::CMASolutions(Parameters &p)
    :_hsig(1),_max_eigenv(0.0),_min_eigenv(0.0),_niter(0),_kcand(1),_eigeniter(0),_updated_eigen(true),_run_status(0),_elapsed_time(0)
  {
    _cov = dMat::Identity(p._dim,p._dim);
    if (p._x0 == -DBL_MAX)
      _xmean = dVec::Random(p._dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
    else _xmean = dVec::Constant(p._dim,p._x0);
    if (static_cast<CMAParameters&>(p)._sigma_init > 0.0)
      _sigma = static_cast<CMAParameters&>(p)._sigma_init;
    else static_cast<CMAParameters&>(p)._sigma_init = _sigma = 1.0/static_cast<double>(p._dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(p._dim);
    _pc = dVec::Zero(p._dim);
    _candidates.resize(p._lambda);
    _kcand = static_cast<int>(1.0+floor(0.1+p._lambda/4.0));
  }

  CMASolutions::~CMASolutions()
  {
  }

  void CMASolutions::update_best_candidates()
  {
    _best_candidates_hist.push_back(_candidates.at(0)); // supposed candidates is sorted.
    _k_best_candidates_hist.push_back(_candidates.at(_kcand));
    
    //debug
    /*std::cerr << "ordered candidates:\n";
    for (size_t i=0;i<_candidates.size();i++)
      {
	std::cerr << _candidates.at(i)._fvalue << " / " << _candidates.at(i)._x.transpose() << std::endl;
	}*/
    //debug
  }

  void CMASolutions::update_eigenv_bounds(const dVec &eigenv)
  {
    _max_eigenv = eigenv.maxCoeff();
    _min_eigenv = eigenv.minCoeff();
  }
  
}
