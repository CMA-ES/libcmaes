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

#ifndef PWQ_BOUND_STRATEGY_H
#define PWQ_BOUND_STRATEGY_H

#include "eo_matrix.h"
#include <vector>

namespace libcmaes
{
  class pwqBoundStrategy
  {
  public:
    pwqBoundStrategy(); // dummy constructor, required for non-pointer default object in GenoPheno.
    pwqBoundStrategy(const double *lbounds, const double *ubounds, const int &dim);
    pwqBoundStrategy(const double *lbounds, const double *ubounds,
		     const double *plbounds, const double *pubounds, const int &dim);
    ~pwqBoundStrategy();

    void to_f_representation(const dVec &x,
			     dVec &y) const;
    
    void to_internal_representation(dVec &x,
				    const dVec &y) const;

    void shift_into_feasible(const dVec &x, dVec &x_s) const;

    double getLBound(const int &k) const { return _lbounds[k]; }
    double getUBound(const int &k) const { return _ubounds[k]; }
    double getPhenoLBound(const int &k) const { return _phenolbounds[k]; }
    double getPhenoUBound(const int &k) const { return _phenoubounds[k]; }

    double getAL(const int &k) const { return _al[k]; }
    double getAU(const int &k) const { return _au[k]; }

    void remove_dimensions(const std::vector<int> &k);

    bool is_id() const
    {
      return _id;
    }
    
  private:
    dVec _lbounds;
    dVec _ubounds;
    dVec _al;
    dVec _au;
    dVec _xlow;
    dVec _xup;
    dVec _r;
    dVec _phenolbounds; /**< differ from _lbounds when another geno/pheno transform applies before bounds. */
    dVec _phenoubounds; /**< differ from _ubounds when another geno/pheno transform applies before bounds. */
    bool _id = false;
  };
}

#endif
