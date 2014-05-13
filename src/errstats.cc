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

#include "errstats.h"
#include <iostream>

namespace libcmaes
{
  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_pk(FitFunc &func,
						 CMAParameters<TGenoPheno> &parameters,
						 const CMASolutions &cmasol,
						 const int &k,
						 const double &vk)
  {
    // copy and re-arrange existing solution.
    CMASolutions ncmasol = cmasol;
    
    // set re-arranged solution as starting point of the new optimization problem.
    ncmasol.reset_as_fixed(k);

    // re-arrange function
    FitFunc nfunc = [func,k,vk](const double *x, const int N)
      {
	std::vector<double> nx;
	nx.reserve(N+1);
	for (int i=0;i<N+1;i++)
	  {
	    if (i < k)
	      nx.push_back(x[i]);
	    else if (i == k)
	      nx.push_back(vk);
	    else nx.push_back(x[i-1]);
	  }
	return func(&nx.front(),N+1);
      };
    
    // optimize and return result.
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters.reset_as_fixed(k);
    return cmaes(nfunc,nparameters);
  }

  template class errstats<GenoPheno<NoBoundStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy>>;
}
