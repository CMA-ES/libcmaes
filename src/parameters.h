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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "eo_matrix.h"
#include "genopheno.h"
#include <string>
#include <ctime>
#include <cmath>
#include <limits>
#include <unordered_map>

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
     * @param x0 initial search point
     * @param lambda number of offsprings sampled at each step
     * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
     * @param gp genotype / phenotype object
     */
  Parameters(const int &dim, const double *x0, const int &lambda=-1,
	     const uint64_t &seed=0, const TGenoPheno &gp=GenoPheno<NoBoundStrategy>())
  :_dim(dim),_lambda(lambda),_seed(seed),_gp(gp)
  {
    if (_lambda == -1) // lambda is unspecified
      _lambda = 4 + floor(3.0*log(_dim));
    if (_seed == 0) // seed is not forced.
      _seed = static_cast<uint64_t>(time(nullptr));
    set_x0(x0);
  }
  
  ~Parameters()
  {
  }

  /**
   * \brief sets initial objective function parameter values to x0 across all dimensions
   * @param x0 initial value
   */
  void set_x0(const double &x0)
  {
    _x0min = _x0max = dVec::Constant(_dim,x0);
  }

  /**
   * \brief sets initial objective function parameter values to array x0
   * @param x0 array of initial parameter values
   */
  void set_x0(const double *x0)
  {
    _x0min = _x0max = dVec(_dim);
    for (int i=0;i<_dim;i++)
      _x0min(i) = _x0max(i) = x0[i];
  }

  /**
   * \brief sets bounds on initial objective function parameter values.
   *        Bounds are the same across all dimensions, and initial value is
   *        sampled uniformly within these bounds.
   * @param x0min lower bound
   * @param x0max upper bound
   */
  void set_x0(const double &x0min, const double &x0max)
  {
    _x0min = dVec::Constant(_dim,x0min);
    _x0max = dVec::Constant(_dim,x0max);
  }

  void set_x0(const dVec &x0)
  {
    _x0min = x0;
    _x0max = x0;
  }
  
  /**
   * \brief sets bounds on initial objective function parameter values.
   *        Initial value is sampled uniformly within these bounds.
   * @param x0min vector of initial lower bounds.
   * @param x0max vector of initial upper bounds.
   */
  void set_x0(const double *x0min, const double *x0max)
  {
    _x0min = dVec(_dim);
    _x0max = dVec(_dim);
    for (int i=0;i<_dim;i++)
      {
	_x0min(i) = x0min[i];
	_x0max(i) = x0max[i];
      }
  }
  
  /**
   * \brief sets bounds on initial objective function parameter values.
   *        Initial value is sampled uniformly within these bounds.
   * @param x0min vector of initial lower bounds.
   * @param x0max vector of initial upper bounds.
   */
  void set_x0(const dVec &x0min, const dVec &x0max)
  {
    _x0min = x0min;
    _x0max = x0max;
  }

  /**
   * \brief freezes a parameter to a given value during optimization.
   * @param index dimension index of the parameter to be frozen
   * @param value frozen value of the parameter
   */
  void set_fixed_p(const int &index, const double &value)
  {
    _fixed_p.insert(std::pair<int,double>(index,value));
  }

  /**
   * \brief sets the maximum number of iterations allowed for the optimization.
   * @param maxiter maximum number of allowed iterations
   */
  void set_max_iter(const int &maxiter)
  {
    _max_iter = maxiter;
  }

  /**
   * \brief sets the maximum budget of function calls allowed for the optimization.
   * @param fevals maximum number of function evaluations
   */
  void set_max_fevals(const int &fevals)
  {
    _max_fevals = fevals;
  }

  /**
   * \brief sets the objective function target value when known.
   * @param val objective function target value
   */
  void set_ftarget(const double &val)
  {
    _ftarget = val;
  }

  /**
   * \brief resets the objective function target value to its inactive state.
   */
  void reset_ftarget()
  {
    _ftarget = -std::numeric_limits<double>::infinity();
  }
  
  int _dim; /**< function space dimensions. */
  int _lambda = -1; /**< number of offsprings. */
  int _max_iter = -1; /**< max iterations. */
  int _max_fevals = -1; /**< max budget as number of function evaluations. */
  
  bool _quiet = false; /**< quiet all outputs. */
  std::string _fplot = ""; /**< plotting file, if specified. */
  dVec _x0min; /**< initial mean vector min bound value for all components. */
  dVec _x0max; /**< initial mean vector max bound value for all components. */
  double _ftarget = -std::numeric_limits<double>::infinity(); /**< optional objective function target value. */
  double _ftargettol = 1e-12; /**< objective function target tolerance. */
  
  uint64_t _seed = 0; /**< seed for random generator. */
  int _algo = 0; /**< selected algorithm. */

  std::unordered_map<int,double> _fixed_p; /**< fixed parameters and values. */
  
  TGenoPheno _gp; /**< genotype / phenotype object. */
  };
  
}

#endif
