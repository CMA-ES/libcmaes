
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

  /**
   * \brief CMA-ES termination criteria, see reference paper in cmastrategy.h
   */
  class CMAStopCriteria
  {
  public:
    /**
     * \brief Constructor: instanciates a predefined set of termination criteria
     *        tests, see reference paper in cmastrategy.h
     */
    CMAStopCriteria();
    ~CMAStopCriteria();

    /**
     * \brief Termination criteria evaluation: the function iterates and 
     *        evaluates the predefined criteria.
     * @return 0 if no termination criteria triggers, the termination code 
     *           otherwise (< 0 for an error, > 1 for a partial success).
     */
    int stop(const CMAParameters &cmap, const CMASolutions &cmas) const;
    
    std::map<int,StopCriteriaFunc> _scriteria; /**< the set of predefined criteria, with priorities. */
    bool _active; /**< whether these termination criteria are active. */
  };
  
}

#endif
