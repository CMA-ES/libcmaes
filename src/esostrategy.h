/**
 * CMA-ES, Covariance Matrix Evolution Strategy
 * Copyright (c) 2014 INRIA
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

#ifndef ESOSTRATEGY_H
#define ESOSTRATEGY_H

#include "eo_matrix.h" // to include Eigen everywhere.
#include "candidate.h"

namespace libcmaes
{
  typedef std::function<double (const double*, const int &n)> FitFunc;

  template<class TParameters,class TSolutions>
    using ProgressFunc = std::function<int (const TParameters, const TSolutions)>; // template aliasing.

  /**
   * \brief Main class describing an evolutionary optimization strategy.
   *        Every algorithm in libcmaes descends from this class, and bring
   *        its functionalities to an ESOptimizer object.
   */
  template<class TParameters,class TSolutions,class TStopCriteria>
    class ESOStrategy
  {
  public:
    /**
     * \brief constructor
     * @param func function to minimize
     * @param parameters optimization parameters
     */
    ESOStrategy(FitFunc &func,
		TParameters &parameters);

  protected:
    ~ESOStrategy();

  public:
    /**
     * \brief Generates a set of candidate points.
     * @return A matrix whose rows contain the candidate points.
     */
    dMat ask();

    /**
     * \brief Evaluates a set of candiates against the objective function.
     *        The procedure is multithreaded and stores both the candidates
     *        and their f-value into the _solutions object that bears the 
     *        current set of potential solutions to the optimization problem.
     * @param candidates A matrix whose rows contain the candidates.
     */
    void eval(const dMat &candidates);

    /**
     * \brief Updates the state of the stochastic search, and prepares
     *        for the next iteration.
     */
    void tell();

    /**
     * \brief Decides whether to stop the search for solutions.
     * @return true if search must stop, false otherwise.
     */
    bool stop();

    /**
     * \brief Finds the minimum of the objective function. It makes
     *        alternative calls to ask(), tell() and stop() until 
     *        one of the termination criteria triggers.
     * @return success or error code, as defined in opti_err.h
     */
    int optimize();

    // deprecated.
    Candidate best_solution() const;
    
    FitFunc _func; /**< the objective function. */
    int _nevals;  /**< number of function evaluations. */
    int _niter;  /**< number of iterations. */
    TSolutions _solutions; /**< holder of the current set of solutions and the dynamic elemenst of the search state in general. */
    TParameters _parameters; /**< the optimizer's set of static parameters, from inputs or internal. */
    ProgressFunc<TParameters,TSolutions> _pfunc; /**< possibly custom progress function. */
  };
  
}

#endif
