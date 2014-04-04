
#include "cmasolutions.h"
#include <iostream>

namespace libcmaes
{

  CMASolutions::CMASolutions(const int &dim,
		       const int &ncandidates)
    :_hsig(1)
  {
    _cov = dMat::Identity(dim,dim);
    _xmean = dVec::Random(dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
    //_xmean = dVec::Constant(dim,0.5);
    _sigma = 0.3;//1.0/static_cast<double>(dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(dim);
    _pc = dVec::Zero(dim);
    _candidates.resize(ncandidates);
  }

  CMASolutions::~CMASolutions()
  {
  }

  void CMASolutions::update_best_candidates()
  {
    _best_candidates_hist.push_back(_candidates.at(0)); // supposed candidates is sorted.

    //debug
    /*std::cerr << "ordered candidates:\n";
    for (size_t i=0;i<_candidates.size();i++)
      {
	std::cerr << _candidates.at(i)._fvalue << " / " << _candidates.at(i)._x.transpose() << std::endl;
	}*/
    //debug
  }
  
}
