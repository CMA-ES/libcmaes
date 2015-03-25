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

#include "libcmaes_config.h"
#include "candidate.h"
#include "eo_matrix.h"
#include "cmaparameters.h"
#include "cmastopcriteria.h"
#include "pli.h"
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
    template <class U, class V> friend class CMAStrategy;
    template <class U, class V, class W> friend class ESOptimizer;
    template <class U, class V, class W> friend class ESOStrategy;
    template <class U> friend class CMAStopCriteria;
    template <class U, class V> friend class IPOPCMAStrategy;
    template <class U, class V> friend class BIPOPCMAStrategy;
    friend class CovarianceUpdate;
    friend class ACovarianceUpdate;
    template <class U> friend class errstats;
#ifdef HAVE_SURROG
    template <template <class X,class Y> class U, class V, class W> friend class SimpleSurrogateStrategy;
    template <template <class X,class Y> class U, class V, class W> friend class ACMSurrogateStrategy;
#endif
    friend class VDCMAUpdate;
    
  public:
    /**
     * \brief dummy constructor.
     */
    CMASolutions() {}

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
      std::stable_sort(_candidates.begin(),_candidates.end(),
		[](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() < c2.get_fvalue();});
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
     * @return current best candidate
     * @see CMASolutions::sort_candidates
     */
    inline Candidate best_candidate() const
    {
      if (_best_candidates_hist.empty()) // iter = 0
	{
	  if (_initial_candidate.get_x_size())
	    return _initial_candidate;
	  else return Candidate(std::numeric_limits<double>::quiet_NaN(),_xmean);
	}
      return _best_candidates_hist.back();
    }

    /**
     * \brief returns the best seen candidate.
     * @return best seen candidate
     */
    inline Candidate get_best_seen_candidate() const
    {
      return _best_seen_candidate;
    }

    /**
     * \brief returns the worst seen candidate.
     * @return worst seen candidate
     */
    inline Candidate get_worst_seen_candidate() const
    {
      return _worst_seen_candidate;
    }

    /**
     * \brief get a reference to the r-th candidate in current set
     * @param r candidate position
     */
    inline Candidate& get_candidate(const int &r)
      {
	return _candidates.at(r);
      }

    inline Candidate get_candidate(const int &r) const
    {
      return _candidates.at(r);
    }

    /**
     * \brief get a reference to the full candidate set
     */
    inline std::vector<Candidate>& candidates()
    {
      return _candidates;
    }
    
    /**
     * \brief number of candidate solutions.
     * @return current number of solution candidates.
     */
    inline int size() const
    {
      return _candidates.size();
    }

    /**
     * \brief resets the solution object in order to restart from
     *        the current solution with fresh covariance matrix.
     * Note: experimental.
     */
    void reset();
    
    /**
     * \brief re-arrange solution object such that parameter 'k' is fixed (i.e. removed).
     * @param k index of the parameter to remove.
     */
    void reset_as_fixed(const int &k);

    /**
     * \brief get profile likelihood if previously computed.
     */
    bool get_pli(const int &k, pli &p) const
    {
      std::map<int,pli>::const_iterator mit;
      if ((mit=_pls.find(k))!=_pls.end())
	{
	  p = (*mit).second;
	  return true;
	}
      return false;
    }

    /**
     * \brief return problem dimension.
     * @return problem dimension
     */
    inline int dim() const
    {
      return _xmean.size();
    }
    
    /**
     * \brief returns expected distance to minimum.
     * @return edm
     */
    inline double edm() const
    {
      return _edm;
    }

    /**
     * \brief returns error covariance matrix
     * @return error covariance matrix
     */
    inline dMat cov() const
    {
      return _cov;
    }

    /**
     * \brief returns reference to error covariance matrix
     * @return error covariance matrix
     */
    inline const dMat& cov_ref() const
    {
      return _cov;
    }
    
    /**
     * \brief returns pointer to covariance matrix array
     * @return pointer to covariance matrix array
     */
    inline const double* cov_data() const
    {
      return _cov.data();
    }
    
    /**
     * \brief returns full covariance matrix. Similar to cov() but in case of linear-sized
     *        algorithms like sep and vd, returns the full covariance matrix anyways.
     * @return full size covariance matrix
     */
    dMat full_cov() const;

    /**
     * \brief returns separable covariance diagonal matrix, only applicable to sep-CMA-ES algorithms.
     * @return error covariance diagonal vector
     */
    inline dMat sepcov() const
    {
      return _sepcov;
    }

    /**
     * \brief returns reference to separable covariance diagonal vector, only applicable to sep-CMA-ES algorithms.
     * @return error covariance diagonal vector
     */
    inline const dMat& sepcov_ref() const
    {
      return _sepcov;
    }
    
    /**
     * \brief returns pointer to covariance diagnoal vector
     * @return pointer to covariance diagonal array
     */
    inline const double* sepcov_data() const
    {
      return _sepcov.data();
    }

    /**
     * \brief returns inverse root square of covariance matrix
     * @return square root of error covariance matrix
     */
    inline dMat csqinv() const
    {
      return _csqinv;
    }

    /**
     * \brief returns inverse root square of separable covariance diagonal matrix, only applicable to sep-CMA-ES algorithms.
     * @return square root of error covariance diagonal matrix
     */
    inline dMat sepcsqinv() const
    {
      return _sepcsqinv;
    }

    /**
     * \returns the unscaled standard deviation vector
     * Note: this is only useful to compare amond standard deviations
     * To get the true rescaled estimate of the error, use errors()
     * @param cmaparams parameter object that hold the genotype/phenotype transform
     * @return unscaled standard deviation vector
     */
    template<class TGenoPheno=GenoPheno<NoBoundStrategy>>
      inline dVec stds(const CMAParameters<TGenoPheno> &cmaparams) const
    {
      dVec phen_xmean = cmaparams.get_gp().pheno(_xmean);
      dVec stds;
      if (!cmaparams.is_sep() && !cmaparams.is_vd())
	stds = _cov.diagonal().cwiseSqrt();
      else if (cmaparams.is_sep())
	stds = _sepcov.cwiseSqrt();
      else if (cmaparams.is_vd())
	stds = (dVec::Constant(cmaparams.dim(),1.0)+_v.cwiseProduct(_v)).cwiseSqrt().cwiseProduct(_sepcov);
      dVec phen_xmean_std = cmaparams.get_gp().pheno(static_cast<dVec>(_xmean + stds));
      return (phen_xmean_std - phen_xmean).cwiseAbs();
    }
    
    /**
     * \returns standard deviation vector
     * @param cmaparams parameter object that hold the genotype/phenotype transform
     * @return standard deviation vector
     */
    template<class TGenoPheno=GenoPheno<NoBoundStrategy>>
      inline dVec errors(const CMAParameters<TGenoPheno> &cmaparams) const
      {
	return std::sqrt(_sigma)*stds(cmaparams);
      }

    /**
     * \brief returns correlation matrix
     * @return correlation matrix
     */
    dMat corr() const;

    /**
     * \brief returns correlation between parameter i and j (useful in large-scale settings)
     * @return correlation between parameter i and j
     */
    double corr(const int &i, const int &j) const;

    /**
     * \brief returns current value of step-size sigma
     * @return current step-size
     */
    inline double sigma() const
    {
      return _sigma;
    }

    /**
     * \brief sets new step-size value, use with care
     * @param sigma step-size value
     */
    inline void set_sigma(const double &sigma)
    {
      _sigma = sigma;
    }

    /**
     * \brief returns current distribution's mean in parameter space
     * @return mean
     */
    inline dVec xmean() const
    {
      return _xmean;
    }

    /**
     * \brief sets the current distributions' mean in parameter space
     * @param xmean mean vector
     */
    inline void set_xmean(const dVec &xmean)
    {
      _xmean = xmean;
    }
    
    /**
     * \brief returns current optimization status.
     * @return status
     */
    inline int run_status() const
    {
      return _run_status;
    }

    /**
     * \brief returns current optimization status' message.
     * @return status message
     */
    inline std::string status_msg() const
      {
	return CMAStopCriteria<>::_scriterias[_run_status];
      }

    /**
     * \brief returns currently elapsed time spent on optimization
     * @return time spent on optimization
     */
    inline int elapsed_time() const
    {
      return _elapsed_time;
    }

    /**
     * \brief returns time spent on last iteration
     * @return time spent on last iteration
     */
    inline int elapsed_last_iter() const
    {
      return _elapsed_last_iter;
    }
    
    /**
     * \brief returns current number of iterations
     * @return number of iterations
     */
    inline int niter() const
    {
      return _niter;
    }

    /**
     * \brief returns current budget (number of objective function calls)
     * @return number of objective function calls
     */
    inline int nevals() const
    {
      return _nevals;
    }
    
    /**
     * \brief returns current minimal eigen value
     * @return minimal eigen value
     */
    inline double min_eigenv() const
    {
      return _min_eigenv;
    }

    /**
     * \brief returns current maximal eigen value
     * @return maximal eigen value
     */
    inline double max_eigenv() const
    {
      return _max_eigenv;
    }

    /**
     * \brief returns whether the last update is lazy
     * @return whether the last update is lazy
     */
    inline bool updated_eigen() const
    {
      return _updated_eigen;
    }

    /**
     * \brief returns current number of objective function evaluations
     * @return number of objective function evaluations
     */
    inline int fevals() const
    {
      return _nevals;
    }

    /**
     * \brief returns last computed eigenvalues
     * @return last computed eigenvalues
     */
    inline dVec eigenvalues() const
    {
      return _leigenvalues;
    }
    
    /**
     * \brief returns last computed eigenvectors
     * @return last computed eigenvectors
     */
    inline dMat eigenvectors() const
    {
      return _leigenvectors;
    }
    
    /**
     * \brief print the solution object out.
     * @param out output stream
     * @param verb_level verbosity level: 0 for short, 1 for debug.
     */
    template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
    std::ostream& print(std::ostream &out,
			const int &verb_level=0,
			const TGenoPheno &gp=TGenoPheno()) const;

  private:
    dMat _cov; /**< covariance matrix. */
    dMat _csqinv; /** inverse root square of covariance matrix. */
    dMat _sepcov;
    dMat _sepcsqinv;
    dVec _xmean; /**< distribution mean. */
    dVec _psigma; /**< cumulation for sigma. */
    dVec _pc; /**< cumulation for covariance. */
    short _hsig = 1; /**< 0 or 1. */
    double _sigma; /**< step size. */
    std::vector<Candidate> _candidates; /**< current set of candidate solutions. */
    std::vector<Candidate> _best_candidates_hist; /**< history of best candidate solutions. */
    int _max_hist = -1; /**< max size of the history, keeps memory requirements fixed. */
    
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

    std::map<int,pli> _pls; /**< profile likelihood for parameters it has been computed for. */
    double _edm = 0.0; /**< expected vertical distance to the minimum. */

    Candidate _best_seen_candidate; /**< best seen candidate along the run. */
    int _best_seen_iter;
    Candidate _worst_seen_candidate;
    Candidate _initial_candidate;
    
    dVec _v; /**< complementary vector for use in vdcma. */

    std::vector<RankedCandidate> _candidates_uh; /**< temporary set of candidates used by uncertainty handling scheme. */
    int _lambda_reev; /**< number of reevaluated solutions at current step. */
    double _suh; /**< uncertainty level computed by uncertainty handling procedure. */
    
    double _tpa_s = 0.0;
    int _tpa_p1 = 0;
    int _tpa_p2 = 1;
    dVec _tpa_x1;
    dVec _tpa_x2;
    dVec _xmean_prev; /**< previous step's mean vector. */
  };

  std::ostream& operator<<(std::ostream &out,const CMASolutions &cmas);
  
}

#endif
