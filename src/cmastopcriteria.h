
#ifndef CMASTOPCRITERIA_H
#define CMASTOPCRITERIA_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include <functional>

namespace libcmaes
{

  typedef std::function<int (const CMAParameters &cmap, const CMASolutions &cmas)> StopCriteriaFunc;

  enum CMAStopCritType
  {
    CONT = 0,
    AUTOMAXITER = -10,
    TOLHISTFUN = 1, // convergence
    EQUALFUNVALS = -11,
    TOLX = -12,
    TOLUPSIGMA = -13,
    STAGNATION = -14,
    CONDITIONCOV = -15
  };

  class CMAStopCriteria
  {
  public:
    CMAStopCriteria();
    ~CMAStopCriteria();

    int stop(const CMAParameters &cmap, const CMASolutions &cmas) const;
    
    std::map<int,StopCriteriaFunc> _scriteria;
    bool _active; /**< whether these termination criteria are active. */
  };
  
}

#endif
