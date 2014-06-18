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

#ifndef CMAPARAMETERS_H
#define CMAPARAMETERS_H

#include "parameters.h"
#include "eo_matrix.h"
#include <float.h>

namespace libcmaes
{
  /**
   * \brief Parameters for various flavors of the CMA-ES algorithm.
   */
  template <class TGenoPheno=GenoPheno<NoBoundStrategy> >
  class CMAParameters : public Parameters<TGenoPheno>
  {
  public:
    CMAParameters() {}; //TODO: var init even if this constructor is not supposed to be used for now.

    /**
     * \brief Constructor.
     * @param dim problem dimensions
     * @param x0 initial search point
     * @param sigma initial distribution step size (positive, otherwise automatically set)
     * @param lambda number of offsprings sampled at each step
     * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
     * @param gp genotype / phenotype object
     * @param sep whether to use sep-CMA-ES, using diagonal covariance matrix (modifies covariance default learning rate)
     */
  CMAParameters(const int &dim,
		const double *x0,
		const double &sigma,
		const int &lambda=-1,
		const uint64_t &seed=0,
		const TGenoPheno &gp=GenoPheno<NoBoundStrategy>());
    ~CMAParameters();
  
    void reset_as_fixed(const int &k);

    /**
     * \brief initialize required parameters based on dim, lambda, x0 and sigma.
     */
    void initialize_parameters();
  
    /**
     * \brief adapt parameters for noisy objective function.
     */
    void set_noisy();
  
    /**
     * \brief fix parameters for sep-CMA-ES, using only the diagonal of covariance matrix.
     */
    void set_sep();

    /**
     * \brief turns stopping criteria MaxIter that automatically stops optimization after a 
     *        number of steps on or off.
     * @param b true or false for turning criteria on or off (on is default in constructor).
     */
    void set_automaxiter(const bool &b) { _has_max_iter = b; }

    /**
     * \brief freezes a parameter to a given value during optimization.
     *        Adapts some generic parameters as well.
     * @param index dimension index of the parameter to be frozen
     * @param value frozen value of the parameter
     */
    void set_fixed_p(const int &index, const double &value);

    /**
     * \brief sets the maximum number of restarts (applies to IPOP and BIPOP).
     * @param nrestarts maximum number of restarts
     */
    void set_restarts(const int &nrestarts)
    {
      _nrestarts = nrestarts;
    }
  
    int _mu; /**< number of candidate solutions used to update the distribution parameters. */
    dVec _weights; /**< offsprings weighting scheme. */
    double _csigma; /**< cumulation constant for step size. */
    double _c1; /**< covariance matrix learning rate for the rank one update using pc. */
    double _cmu; /**< covariance matrix learning reate for the rank mu update. */
    double _cc; /**< cumulation constant for pc. */
    double _muw; /**< \sum^\mu _weights .*/
    double _dsigma; /**< step size damping factor. */
    
    // computed once at init for speeding up operations.
    double _fact_ps;
    double _fact_pc;
    double _chi; /**< norm of N(0,I) */

    double _sigma_init; /**< initial sigma value. */

    int _nrestarts = 9; /**< maximum number of restart, when applicable. */
    bool _lazy_update; /**< covariance lazy update. */
    double _lazy_value; /**< reference trigger for lazy update. */

    // active cma.
    double _cm; /**< learning rate for the mean. */
    double _alphacov; /**< = 2 (active CMA only) */
    double _alphaminusold; /**< in [0,1] (active CMA only) */
    double _deltamaxsigma; /**< infinite (active CMA only) */
    double _lambdamintarget; /**< = 0.66 (active CMA only) */
    double _alphaminusmin; /**< = 1 (active CMA only) */

    // sep cma (diagonal cov).
    bool _sep = false; /**< whether to use diagonal covariance matrix. */

    // stopping criteria.
    bool _has_max_iter = true; /**< MaxIter criteria: automatically stop running after 100+50*((D+2)^2)/lambda iterations. */
  };
  
}

#endif
