
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
    AUTOMAXITER = 10,
    TOLHISTFUN = 11,
    EQUALFUNVALS = 12,
    TOLX = 13,
    TOLUPSIGMA = 14,
    STAGNATION = 15,
    CONDITIONCOV = 16
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
