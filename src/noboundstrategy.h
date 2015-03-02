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

#ifndef NOBOUNDSTRATEGY_H
#define NOBOUNDSTRATEGY_H

#include "eo_matrix.h"
#include <vector>
#include <limits>

namespace libcmaes
{
  class NoBoundStrategy
  {
  public:
    NoBoundStrategy(const double *lbounds=nullptr,const double *ubounds=nullptr,const int dim=0)
      {
	(void)lbounds;
	(void)ubounds;
	(void)dim;
      }; // empty constructor with signature.

    NoBoundStrategy(const double *lbounds,const double *ubounds,
		    const double *plbounds,const double *pubounds,const int dim=0)
      {
	(void)lbounds;
	(void)ubounds;
	(void)plbounds;
	(void)pubounds;
	(void)dim;
      }; // empty constructor with signature.
    
    ~NoBoundStrategy() {};

    void to_f_representation(const dVec &x, dVec &y) const
    {
      (void)x;
      (void)y;
    }

    void remove_dimensions(const std::vector<int> &k)
    {
      (void)k;
    }

    bool is_id() const
    {
      return _id;
    }

    double getLBound(const int &k) const { (void)k;return -std::numeric_limits<double>::max(); }
    double getUBound(const int &k) const { (void)k;return std::numeric_limits<double>::max(); }
    double getPhenoLBound(const int &k) const { (void)k;return -std::numeric_limits<double>::max(); }
    double getPhenoUBound(const int &k) const { (void)k;return std::numeric_limits<double>::max(); }
    
  private:
    bool _id = true;
  };
}

#endif
