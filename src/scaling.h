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

/**
 * \brief linear scaling of the parameter space to achieve similar sensitivity
 *        across all components.
 */

#ifndef LSCALING_H
#define LSCALING_H

#include "eo_matrix.h"
#include <limits>
#include <iostream>

namespace libcmaes
{

  class NoScalingStrategy
  {
    friend class CMASolutions;
    template <class U, class V> friend class GenoPheno;
    
  public:
    NoScalingStrategy() {}

    NoScalingStrategy(const double *lbounds,
		      const double *ubounds,
		      const int &dim)
      {
	(void)lbounds;
	(void)ubounds;
	(void)dim;
      }
    
    ~NoScalingStrategy() {}

    void scale_to_internal(dVec &x,
			   const dVec &y) const
    {
      x = y;
    }

    void scale_to_f(const dVec &x,
		    dVec &y) const
    {
      y = x;
    }

    void remove_dimensions(const std::vector<int> &k)
    {
      (void)k;
    }

    bool is_id() const
    {
      return _id;
    }

  private:
    double _intmin = -std::numeric_limits<double>::max();  /**< default internal min bound. */
    double _intmax = std::numeric_limits<double>::max();  /**< default internal max bound. */
    bool _id = true;
  };
  
  class linScalingStrategy
  {
    friend class CMASolutions;
    template <class U, class V> friend class GenoPheno;
    
  public:
    linScalingStrategy() // identity scaling
      :_scaling(dVec::Constant(1,1.0)),_shift(dVec::Zero(1)),_id(true)
      {
      }
    
    linScalingStrategy(const double *lbounds,
		       const double *ubounds,
		       const int &dim)
      :_id(false)
      {
	compute_scaling(lbounds,ubounds,dim);
      }
    
    linScalingStrategy(const dVec &scaling,
		       const dVec &shift)
      :_scaling(scaling),_shift(shift),_id(false)
      {}
    
    ~linScalingStrategy() {}

    void compute_scaling(const double *lbounds,
			 const double *ubounds,
			 const int &dim)
    {
      dVec vlbounds = Eigen::Map<dVec>(const_cast<double*>(lbounds),dim);
      dVec vubounds = Eigen::Map<dVec>(const_cast<double*>(ubounds),dim);
      dVec denom = vubounds-vlbounds;
      denom = denom.cwiseMin(std::numeric_limits<double>::max()); // protects against overflow
      _scaling = (dVec::Constant(dim,_intmax)-dVec::Constant(dim,_intmin)).cwiseQuotient(denom);
      _shift = dVec::Constant(dim,_intmax) - _scaling.cwiseProduct(vubounds);
    }
    
    void scale_to_internal(dVec &x,
			   const dVec &y) const
    {
      x = y.cwiseProduct(_scaling) + _shift;
    }

    void scale_to_f(const dVec &x,
		    dVec &y) const
    {
      y = x - _shift;
      y = y.cwiseQuotient(_scaling);
    }

    bool is_id() const { return _id; }

    void remove_dimensions(const std::vector<int> &k)
    {
      for (const int i: k)
	{
	  removeElement(_scaling,i);
	  removeElement(_shift,i);
	}
    }
    
  public:
    double _intmin = 0.0;  /**< default internal min bound. */
    double _intmax = 10.0;  /**< default internal max bound. */
    dVec _scaling;
    dVec _shift;

  public:
    bool _id = true;
  };
  
}

#endif
