
#ifndef CMAES_H
#define CAMES_H

#include "esoptimizer.h"
#include "cmastrategy.h"
#include "ipopcmastrategy.h"
#include "bipopcmastrategy.h"

namespace libcmaes
{
  CMASolutions cmaes(FitFunc &func,
		     CMAParameters &parameters,
		     ProgressFunc<CMAParameters,CMASolutions> &pfunc=CMAStrategy<CovarianceUpdate>::_defaultPFunc)
  {
    if (parameters._algo == CMAES_DEFAULT)
      {
	ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters> cmaes_vanilla(func,parameters);
	cmaes_vanilla.set_progress_func(pfunc);
	cmaes_vanilla.optimize();
	return cmaes_vanilla._solutions;
      }
    else if (parameters._algo == IPOP_CMAES)
      {
	ESOptimizer<IPOPCMAStrategy,CMAParameters> ipop(func,parameters);
	ipop.set_progress_func(pfunc);
	ipop.optimize();
	return ipop._solutions;
      }
    else if (parameters._algo == BIPOP_CMAES)
      {
	ESOptimizer<BIPOPCMAStrategy,CMAParameters> bipop(func,parameters);
	bipop.set_progress_func(pfunc);
	bipop.optimize();
	return bipop._solutions;
      }
    return CMASolutions();
  };
}

#endif
