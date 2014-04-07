
#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "eo_matrix.h"

namespace libcmaes
{
  
  class Candidate
  {
  public:
  Candidate():
    _fvalue(0.0) {};
  Candidate(const double &fvalue,
	    const dVec &x)
    :_fvalue(fvalue),_x(x)
    {};
  ~Candidate() {};
    
    double _fvalue; /**< function value. */
    dVec _x;
  };

}

#endif
