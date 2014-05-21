/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
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

#ifndef ERRSTATS_H
#define ERRSTATS_H

#include "pli.h"
#include "cmaes.h"

namespace libcmaes
{
  
  template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class errstats
    {
    public:
    static pli profile_likelihood(FitFunc &func,
				  const CMAParameters<TGenoPheno> &parameters,
				  CMASolutions &cmasol,
				  const int &k,
				  const bool &curve=false,
				  const int &samplesize=1000,
				  const double &fup=0.1,
				  const double &delta=0.1);

    static void profile_likelihood_search(FitFunc &func,
					  const CMAParameters<TGenoPheno> &parameters,
					  pli &le,
					  const CMASolutions &cmasol,
					  const int &k,
					  const bool &neg,
					  const int &samplesize,
					  const double &fup,
					  const double &delta,
					  const bool &curve);
    
    static bool take_linear_step(FitFunc &func,
				 const CMAParameters<TGenoPheno> &parameters,
				 const int &k,
				 const double &minfvalue,
				 const double &fup,
				 const bool &curve,
				 dVec &x,
				 double &dxk);

    static CMASolutions optimize_pk(FitFunc &func,
				    const CMAParameters<TGenoPheno> &parameters,
				    const CMASolutions &cmasol,
				    const int &k,
				    const double &vk);
    
    // DEPRECATED
    static CMASolutions optimize_reduced_pk(FitFunc &func,
					    CMAParameters<TGenoPheno> &parameters,
					    const CMASolutions &cmasol,
					    const int &k,
					    const double &vk);
    };    
  
}

#endif
