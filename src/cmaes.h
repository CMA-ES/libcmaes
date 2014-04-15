
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
	ESOptimizer<IPOPCMAStrategy<CovarianceUpdate>,CMAParameters> ipop(func,parameters);
	ipop.set_progress_func(pfunc);
	ipop.optimize();
	return ipop._solutions;
      }
    else if (parameters._algo == BIPOP_CMAES)
      {
	ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate>,CMAParameters> bipop(func,parameters);
	bipop.set_progress_func(pfunc);
	bipop.optimize();
	return bipop._solutions;
      }
    else if (parameters._algo == aCMAES)
      {
	ESOptimizer<CMAStrategy<ACovarianceUpdate>,CMAParameters> acmaes(func,parameters);
	acmaes.set_progress_func(pfunc);
	acmaes.optimize();
	return acmaes._solutions;
      }
    else if (parameters._algo == aIPOP_CMAES)
      {
	ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate>,CMAParameters> aipop(func,parameters);
	aipop.set_progress_func(pfunc);
	aipop.optimize();
	return aipop._solutions;
      }
    else if (parameters._algo == aBIPOP_CMAES)
      {
	ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate>,CMAParameters> abipop(func,parameters);
	abipop.set_progress_func(pfunc);
	abipop.optimize();
	return abipop._solutions;
      }
    return CMASolutions();
  };
}

#endif
