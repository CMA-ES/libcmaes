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

#ifndef CMAES_H
#define CAMES_H

#include "esoptimizer.h"
#include "cmastrategy.h"
#include "ipopcmastrategy.h"
#include "bipopcmastrategy.h"

namespace libcmaes
{
  CMASolutions cmaes(FitFunc &func,
		     CMAParameters<NoBoundStrategy> &parameters,
		     ProgressFunc<CMAParameters<NoBoundStrategy>,CMASolutions> &pfunc=CMAStrategy<CovarianceUpdate>::_defaultPFunc)
  {
    if (parameters._algo == CMAES_DEFAULT)
      {
	ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<NoBoundStrategy>> cmaes_vanilla(func,parameters);
	cmaes_vanilla.set_progress_func(pfunc);
	cmaes_vanilla.optimize();
	return cmaes_vanilla._solutions;
      }
    else if (parameters._algo == IPOP_CMAES)
      {
	ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,NoBoundStrategy>,CMAParameters<NoBoundStrategy>> ipop(func,parameters);
	ipop.set_progress_func(pfunc);
	ipop.optimize();
	return ipop._solutions;
      }
    else if (parameters._algo == BIPOP_CMAES)
      {
	ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,NoBoundStrategy>,CMAParameters<NoBoundStrategy>> bipop(func,parameters);
	bipop.set_progress_func(pfunc);
	bipop.optimize();
	return bipop._solutions;
      }
    else if (parameters._algo == aCMAES)
      {
	ESOptimizer<CMAStrategy<ACovarianceUpdate,NoBoundStrategy>,CMAParameters<NoBoundStrategy>> acmaes(func,parameters);
	acmaes.set_progress_func(pfunc);
	acmaes.optimize();
	return acmaes._solutions;
      }
    else if (parameters._algo == aIPOP_CMAES)
      {
	ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,NoBoundStrategy>,CMAParameters<NoBoundStrategy>> aipop(func,parameters);
	aipop.set_progress_func(pfunc);
	aipop.optimize();
	return aipop._solutions;
      }
    else if (parameters._algo == aBIPOP_CMAES)
      {
	ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,NoBoundStrategy>,CMAParameters<NoBoundStrategy>> abipop(func,parameters);
	abipop.set_progress_func(pfunc);
	abipop.optimize();
	return abipop._solutions;
      }
    return CMASolutions();
  }
}

#endif
