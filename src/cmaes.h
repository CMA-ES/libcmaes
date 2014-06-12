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
		     ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> &pfunc=CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultPFunc,
		     const CMASolutions &solutions=CMASolutions())
    {
      switch(parameters._algo)
	{
	case CMAES_DEFAULT:
	{
	  // if / else as object cannot be reused for now, and we're avoiding pointers.
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> cmaes_vanilla(func,parameters);
	      cmaes_vanilla.set_progress_func(pfunc);
	      cmaes_vanilla.optimize();
	      return cmaes_vanilla._solutions;
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> cmaes_vanilla(func,parameters,solutions);
	      cmaes_vanilla.set_progress_func(pfunc);
	      cmaes_vanilla.optimize();
	      return cmaes_vanilla._solutions;
	    }
	}
	case IPOP_CMAES:
	{
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      ipop.set_progress_func(pfunc);
	      ipop.optimize();
	      return ipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters,solutions);
	      ipop.set_progress_func(pfunc);
	      ipop.optimize();
	      return ipop._solutions;
	    }
	}
	case BIPOP_CMAES:
	{
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters,solutions);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop._solutions;
	    }
	}
	case aCMAES:
	{
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> acmaes(func,parameters);
	      acmaes.set_progress_func(pfunc);
	      acmaes.optimize();
	      return acmaes._solutions;
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> acmaes(func,parameters,solutions);
	      acmaes.set_progress_func(pfunc);
	      acmaes.optimize();
	      return acmaes._solutions;
	    }
	}
	case aIPOP_CMAES:
	{
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> aipop(func,parameters);
	      aipop.set_progress_func(pfunc);
	      aipop.optimize();
	      return aipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> aipop(func,parameters,solutions);
	      aipop.set_progress_func(pfunc);
	      aipop.optimize();
	      return aipop._solutions;
	    }
	}
	case aBIPOP_CMAES:
	{
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> abipop(func,parameters);
	      abipop.set_progress_func(pfunc);
	      abipop.optimize();
	      return abipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> abipop(func,parameters,solutions);
	      abipop.set_progress_func(pfunc);
	      abipop.optimize();
	      return abipop._solutions;
	    }
	}
	case sepCMAES:
	{
	  parameters.set_sep();
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.optimize();
	      return sepcmaes._solutions;
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.optimize();
	      return sepcmaes._solutions;
	    }
	}
	case sepIPOP_CMAES:
	{
	  parameters.set_sep();
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      ipop.set_progress_func(pfunc);
	      ipop.optimize();
	      return ipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      ipop.set_progress_func(pfunc);
	      ipop.optimize();
	      return ipop._solutions;
	    }
	}
	case sepBIPOP_CMAES:
	{
	  parameters.set_sep();
	  if (solutions._cov.size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop._solutions;
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop._solutions;
	    }
	}
	case sepaCMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters);
	  sepcmaes.set_progress_func(pfunc);
	  sepcmaes.optimize();
	  return sepcmaes._solutions;
	}
	case sepaIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	  ipop.set_progress_func(pfunc);
	  ipop.optimize();
	  return ipop._solutions;
	}
	case sepaBIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	  bipop.set_progress_func(pfunc);
	  bipop.optimize();
	  return bipop._solutions;
	}
	default:
	return CMASolutions();
	}
    }
}

#endif
