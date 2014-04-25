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

#ifndef PWQ_BOUND_STRATEGY_H
#define PWQ_BOUND_STRATEGY_H

#include "eo_matrix.h"

//TODO: NoBounds class.
//TODO: make a lambda function typedef for the mapping functions (should be in genopheno...)
//TODO: unit tests.

namespace libcmaes
{
  class pwqBoundStrategy
  {
  public:
    pwqBoundStrategy(const double *lbounds, const double *ubounds, const int &dim);
    ~pwqBoundStrategy();

    void to_f_representation(const dVec &x,
			     dVec &y);
    
    void to_internal_representation(dVec &x,
				    const dVec &y);

    void shift_into_feasible(const dVec &x, dVec &x_s);
    
  public:
    dVec _lbounds;
    dVec _ubounds;
    dVec _al;
    dVec _au;
    dVec _xlow;
    dVec _xup;
    dVec _r;
  };
}

#endif
