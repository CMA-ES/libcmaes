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

#ifndef SURRCMAES_H
#define SURRCMAES_H

#include "cmaes.h"
#include "surrogates/rankingsvm.hpp"
#include "surrogates/rsvm_surr_strategy.hpp"

namespace libcmaes
{
  template <class TGenoPheno=GenoPheno<NoBoundStrategy> >
  CMASolutions CMAES_EXPORT surrcmaes(FitFunc &func,
			 CMAParameters<TGenoPheno> &parameters)
    {
      switch(parameters.get_algo())
	{
	case CMAES_DEFAULT:
	{
	  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case IPOP_CMAES:
	{
	  ESOptimizer<RSVMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case BIPOP_CMAES:
	{
	  ESOptimizer<RSVMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case aCMAES:
	{
	  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case aIPOP_CMAES:
	{
	  ESOptimizer<RSVMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case aBIPOP_CMAES:
	{
	  ESOptimizer<RSVMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepCMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepBIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepaCMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepaIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case sepaBIPOP_CMAES:
	{
	  parameters.set_sep();
	  ESOptimizer<RSVMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case VD_CMAES:
	{
	  parameters.set_vd();
	  ESOptimizer<RSVMSurrogateStrategy<CMAStrategy,VDCMAUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case VD_IPOP_CMAES:
	{
	  parameters.set_vd();
	  ESOptimizer<RSVMSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	case VD_BIPOP_CMAES:
	{
	  parameters.set_vd();
	  ESOptimizer<RSVMSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate>,CMAParameters<>> optim(func,parameters);
	  optim.optimize();
	  return optim.get_solutions();
	}
	default:
	return CMASolutions();
	}
    }
}

#endif
