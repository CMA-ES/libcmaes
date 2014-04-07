
#include "cmasolutions.h"
#include <iostream>

namespace libcmaes
{

  CMASolutions::CMASolutions(const int &dim,
		       const int &ncandidates)
    :_hsig(1),_max_eigenv(0.0),_min_eigenv(0.0),_niter(0),_kcand(1)
  {
    _cov = dMat::Identity(dim,dim);
    _xmean = dVec::Random(dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
    //_xmean = dVec::Constant(dim,0.5);
    _sigma = _sigma_init = 1.0/static_cast<double>(dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(dim);
    _pc = dVec::Zero(dim);
    _candidates.resize(ncandidates);
    _kcand = static_cast<int>(1.0+floor(0.1+ncandidates/4.0));
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
