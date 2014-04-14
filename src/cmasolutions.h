
#ifndef CMASOLUTIONS_H
#define CMASOLUTIONS_H

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
     * DO NOT USE.
     */
    CMASolutions() {};

    /**
     * \brief initializes solutions from stochastic optimization parameters.
     * @param p parameters
     */
    CMASolutions(Parameters &p);
    
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
    short _hsig; /**< 0 or 1. */
    double _sigma; /**< step size. */
    std::vector<Candidate> _candidates; /**< current set of candidate solutions. */
    std::vector<Candidate> _best_candidates_hist; /**< history of best candidate solutions. */

    double _max_eigenv; /**< max eigenvalue, for termination criteria. */
    double _min_eigenv; /**< min eigenvalue, for termination criteria. */
    dVec _leigenvalues; /**< last computed eigenvalues, for termination criteria. */
    dMat _leigenvectors; /**< last computed eigenvectors, for termination criteria. */
    int _niter; /**< number of iterations to reach this solution, for termination criteria. */
    int _kcand;
    std::vector<Candidate> _k_best_candidates_hist; /**< k-th best candidate history, for termination criteria, k is kcand=1+floor(0.1+lambda/4). */
    std::vector<double> _bfvalues; /**< best function values over the past 20 steps, for termination criteria. */
    std::vector<double> _median_fvalues; /**< median function values of some steps, in the past, for termination criteria. */
    
    int _eigeniter; /**< eigenvalues computation last step, lazy-update only. */
    bool _updated_eigen; /**< last update is not lazy. */

    // status of the run.
    int _run_status; /**< current status of the stochastic optimization (e.g. running, or stopped under termination criteria). */
    int _elapsed_time; /**< final elapsed time of stochastic optimization. */
  };

  std::ostream& operator<<(std::ostream &out,const CMASolutions &cmas);
  
}

#endif
