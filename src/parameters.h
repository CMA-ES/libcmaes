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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "eo_matrix.h"
#include "genopheno.h"
#include <string>
#include <time.h>
#include <math.h>
#include <limits>

namespace libcmaes
{
  /**
   * \brief Generic class for Evolution Strategy parameters.
   */
  template <class TGenoPheno=GenoPheno<NoBoundStrategy> >
  class Parameters
  {
  public:
    /**
     * \brief empty constructor.
     */
  Parameters():_dim(0),_lambda(0),_max_iter(0)
      {}

    /**
     * \brief constructor
     * @param dim problem dimensions
     * @param lambda number of offsprings sampled at each step
     * @param max_iter maximum number of iterations
     * @param max_fevals function evaluation budget as the max number of calls
     * @param x0 initial value of the search in parameter space (if unspecified, sampled from within [-4,4] in all coordinates)
     * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
     */
  Parameters(const int &dim, const int &lambda=-1, const int &max_iter=-1,
	     const int &max_fevals=-1,
	     const double &x0=std::numeric_limits<double>::min(), const std::string &fplot="",
	     const uint64_t &seed=0)
    :_dim(dim),_lambda(lambda),_max_iter(max_iter),_max_fevals(max_fevals),
      _quiet(false),_fplot(fplot),_seed(seed),_algo(0)
      {
	if (_lambda == -1) // lambda is unspecified
	  _lambda = 4 + floor(3.0*log(_dim));
      if (_seed == 0) // seed is not forced.
	_seed = static_cast<uint64_t>(time(NULL));
      _x0min = _x0max = dVec::Constant(_dim,x0);
    }
  Parameters(const int &dim, const int &lambda, const int &max_iter,
	     const int &max_fevals=-1,
	     const double &x0min=std::numeric_limits<double>::min(),
	     const double &x0max=std::numeric_limits<double>::max(),
	     const std::string &fplot="",
	     const uint64_t &seed=0)
    :_dim(dim),_lambda(lambda),_max_iter(max_iter),_max_fevals(max_fevals),
      _quiet(false),_fplot(fplot),_x0min(x0min),_x0max(x0max),_seed(seed),_algo(0)
    {
      if (_seed == 0) // seed is not forced.
	_seed = static_cast<uint64_t>(time(NULL));
    }
  ~Parameters()
    {
    }

    int _dim; /**< function space dimensions. */
    int _lambda; /**< number of offsprings. */
    int _max_iter; /**< max iterations. */
    int _max_fevals; /**< max budget as number of function evaluations. */
    
    bool _quiet; /**< quiet all outputs. */
    std::string _fplot; /**< plotting file, if specified. */
    dVec _x0min; /**< initial mean vector min bound value for all components. */
    dVec _x0max; /**< initial mean vector max bound value for all components. */
    
    uint64_t _seed; /**< seed for random generator. */
    int _algo; /**< selected algorithm. */

    TGenoPheno _gp;
  };
  
}

#endif
