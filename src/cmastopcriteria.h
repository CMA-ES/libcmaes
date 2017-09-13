/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CMASTOPCRITERIA_H
#define CMASTOPCRITERIA_H

#include "cmaparameters.h"
#include <functional>
#include <map>

namespace libcmaes
{

  class CMASolutions;

  template <class TGenoPheno>
    using StopCriteriaFunc = std::function<int (const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)>;

  enum CMAStopCritType
  {
    CONT = 0,
    AUTOMAXITER = 7,
    TOLHISTFUN = 1, // convergence
    EQUALFUNVALS = 5, // partial success, user error
    TOLX = 2, // partial success
    TOLUPSIGMA = -13,
    STAGNATION = 6, // partial success
    CONDITIONCOV = -15, // error, user action needed
    NOEFFECTAXIS = 3, // partial success
    NOEFFECTCOOR = 4, // partial success
    MAXFEVALS = 8,
    MAXITER = 9,
    FTARGET = 10 // success
  };

  template <class TGenoPheno=NoBoundStrategy>
  class StopCriteria
  {
  public:
    StopCriteria(const StopCriteriaFunc<TGenoPheno> &sfunc)
      :_sfunc(sfunc)
    {}
    ~StopCriteria() {}

    inline bool active() const { return _active; }
  
    void set_active(const bool &a) { _active = a; }
  
    StopCriteriaFunc<TGenoPheno> _sfunc;
    bool _active = true;
  };
  
  /**
   * \brief CMA-ES termination criteria, see reference paper in cmastrategy.h
   */
  template <class TGenoPheno=NoBoundStrategy>
  class CMAES_EXPORT CMAStopCriteria
  {
    friend class CMASolutions;

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
    int stop(const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas) const;

    /**
     * \brief activates / deactivates a stopping criteria
     * @param c the criteria to modify
     * @param true to activate, false to deactivate
     * @return 1 if criteria cannot be found, 0 otherwise
     */
    int set_criteria_active(const int &c, const bool &active);
  
  private:
    std::map<int,StopCriteria<TGenoPheno> > _scriteria; /**< the set of predefined criteria, with priorities. */
    bool _active; /**< whether these termination criteria are active. */
    static std::map<int,std::string> _scriterias;
  };

}

#endif
