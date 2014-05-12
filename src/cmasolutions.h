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

#ifndef CMASOLUTIONS_H
#define CMASOLUTIONS_H

#include "config.h"
#include "candidate.h"
#include "eo_matrix.h"
#include "cmaparameters.h"
#include <vector>
#include <algorithm>

namespace libcmaes
{
  /**
   * \brief Holder of the set of evolving solutions from running an instance
   *        of CMA-ES.
   */
  class CMASolutions
  {
  public:
    /**
     * \brief dummy constructor.
     */
    CMASolutions() {};

    /**
     * \brief initializes solutions from stochastic optimization parameters.
     * @param p parameters
     */
    template<class TGenoPheno=GenoPheno<NoBoundStrategy>>
    CMASolutions(Parameters<TGenoPheno> &p);
    
    ~CMASolutions();

    /**
     * \brief sorts the current internal set of solution candidates.
     */
    void sort_candidates()
    {
      std::sort(_candidates.begin(),_candidates.end(),
		[](Candidate const &c1, Candidate const &c2){return c1._fvalue < c2._fvalue;});
    }

    /**
     * \brief updates the history of best candidates, as well as other meaningful
     *        values, typically used in termination criteria.
     * @see CMAStopCriteria
     */
    void update_best_candidates();

    /**
     * \brief updates reference eigenvalue and eigenvectors, for use in 
     *        termination criteria.
     * @see CMAStopCriteria
     */
    void update_eigenv(const dVec &eigenvalues,
		       const dMat &eigenvectors);

    /**
     * \brief returns current best solution candidate.
     *        NOTE: candidates MUST be sorted
     * @return currentbest candidate
     * @see CMASolutions::sort_candidates
     */
    Candidate best_candidate() const
    {
      return _best_candidates_hist.back();
    }

    /**
     * \brief number of candidate solutions.
     * @return current number of solution candidates.
     */
    int size() const
    {
      return _candidates.size();
    }

    /**
     * \brief print the solution object out.
     * @param out output stream
     * @param verb_level verbosity level: 0 for short, 1 for debug.
     */
    std::ostream& print(std::ostream &out,
			const int &verb_level=0) const;

    dMat _cov; /**< covariance matrix. */
    dMat _csqinv; /** inverse root square of covariance matrix. */
    dVec _xmean; /**< distribution mean. */
    dVec _psigma; /**< cummulation for sigma. */
    dVec _pc; /**< cumulation for covariance. */
    short _hsig = 1; /**< 0 or 1. */
    double _sigma; /**< step size. */
    std::vector<Candidate> _candidates; /**< current set of candidate solutions. */
    std::vector<Candidate> _best_candidates_hist; /**< history of best candidate solutions. */

    double _max_eigenv = 0.0; /**< max eigenvalue, for termination criteria. */
    double _min_eigenv = 0.0; /**< min eigenvalue, for termination criteria. */
    dVec _leigenvalues; /**< last computed eigenvalues, for termination criteria. */
    dMat _leigenvectors; /**< last computed eigenvectors, for termination criteria. */
    int _niter = 0; /**< number of iterations to reach this solution, for termination criteria. */
    int _nevals = 0; /**< number of function calls to reach the current solution. */
    int _kcand = 1;
    std::vector<Candidate> _k_best_candidates_hist; /**< k-th best candidate history, for termination criteria, k is kcand=1+floor(0.1+lambda/4). */
    std::vector<double> _bfvalues; /**< best function values over the past 20 steps, for termination criteria. */
    std::vector<double> _median_fvalues; /**< median function values of some steps, in the past, for termination criteria. */
    
    int _eigeniter = 0; /**< eigenvalues computation last step, lazy-update only. */
    bool _updated_eigen = true; /**< last update is not lazy. */

    // status of the run.
    int _run_status = 0; /**< current status of the stochastic optimization (e.g. running, or stopped under termination criteria). */
    int _elapsed_time = 0; /**< final elapsed time of stochastic optimization. */
    int _elapsed_last_iter = 0; /**< time consumed during last iteration. */
#ifdef HAVE_DEBUG
    int _elapsed_eval = 0;
    int _elapsed_ask = 0;
    int _elapsed_tell = 0;
    int _elapsed_stop = 0;
#endif
  };

  std::ostream& operator<<(std::ostream &out,const CMASolutions &cmas);
  
}

#endif
