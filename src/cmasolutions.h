
#ifndef CMASOLUTIONS_H
#define CMASOLUTIONS_H

#include "candidate.h"
#include "eo_matrix.h"
#include <vector>
#include <algorithm>

namespace libcmaes
{
  class CMASolutions
  {
  public:
    CMASolutions() {};
    CMASolutions(const int &dim,
		 const int &ncandidates);
    ~CMASolutions();

    void sort()
    {
      std::sort(_candidates.begin(),_candidates.end(),
		[](Candidate const &c1, Candidate const &c2){return c1._fvalue < c2._fvalue;});
    }

    void update_best_candidates();

    void update_eigenv_bounds(const dVec &eigenv);
    
    Candidate best_candidate() const
    {
      return _best_candidates_hist.back();
    }
    
    int size() const
    {
      return _candidates.size();
    }
    
    dMat _cov; /**< covariance matrix. */
    dVec _xmean; /**< distribution mean. */
    dVec _psigma;
    dVec _pc;
    short _hsig; /**< 0 or 1. */
    double _sigma;
    std::vector<Candidate> _candidates;
    std::vector<Candidate> _best_candidates_hist;

    double _sigma_init; /**< initial sigma, used in termination criteria. */
    double _max_eigenv; /**< max eigenvalue, for termination criteria. */
    double _min_eigenv; /**< min eigenvalue, for termination criteria. */
    int _niter; /**< number of iterations to reach this solution, for termination criteria. */
    int _kcand;
    std::vector<Candidate> _k_best_candidates_hist; /**< k-th best candidate history, for termination criteria, k is kcand=1+floor(0.1+lambda/4). */
  };

}

#endif
