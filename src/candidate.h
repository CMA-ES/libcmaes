
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
    ~Candidate() {};
    
    double _fvalue; /**< function value. */
    dVec _x;
  };

}

#endif
