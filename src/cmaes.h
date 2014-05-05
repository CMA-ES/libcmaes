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

#ifndef CMAES_H
#define CMAES_H

#include "esoptimizer.h"
#include "cmastrategy.h"
#include "ipopcmastrategy.h"
#include "bipopcmastrategy.h"

namespace libcmaes
{
  template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  CMASolutions cmaes(FitFunc &func,
		     CMAParameters<TGenoPheno> &parameters,
		     ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> &pfunc=CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultPFunc)
    {
      switch(parameters._algo)
	{
	case CMAES_DEFAULT:
	{
	  ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> cmaes_vanilla(func,parameters);
	  cmaes_vanilla.set_progress_func(pfunc);
	  cmaes_vanilla.optimize();
	  return cmaes_vanilla._solutions;
	}
	case IPOP_CMAES:
	{
	  ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> ipop(func,parameters);
	  ipop.set_progress_func(pfunc);
	  ipop.optimize();
	  return ipop._solutions;
	}
	case BIPOP_CMAES:
	{
	  ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> bipop(func,parameters);
	  bipop.set_progress_func(pfunc);
	  bipop.optimize();
	  return bipop._solutions;
	}
	case aCMAES:
	{
	  ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> acmaes(func,parameters);
	  acmaes.set_progress_func(pfunc);
	  acmaes.optimize();
	  return acmaes._solutions;
	}
	case aIPOP_CMAES:
	{
	  ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> aipop(func,parameters);
	  aipop.set_progress_func(pfunc);
	  aipop.optimize();
	  return aipop._solutions;
	}
	case aBIPOP_CMAES:
	{
	  ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> abipop(func,parameters);
	  abipop.set_progress_func(pfunc);
	  abipop.optimize();
	  return abipop._solutions;
	}
	default:
	return CMASolutions();
	}
    }
}

#endif
