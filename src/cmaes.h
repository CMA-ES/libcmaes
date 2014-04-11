
#ifndef CMAES_H
#define CAMES_H

#include "esoptimizer.h"
#include "cmastrategy.h"
#include "ipopcmastrategy.h"

namespace libcmaes
{
  //TODO: return solution object.
  CMASolutions cmaes(FitFunc &func,
		     CMAParameters &parameters)
  {
    if (parameters._algo == CMAES_DEFAULT)
      {
	ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters> cmaes_vanilla(func,parameters);
	cmaes_vanilla.optimize();
	return cmaes_vanilla._solutions;
      }
    else if (parameters._algo == IPOP_CMAES)
      {
	ESOptimizer<IPOPCMAStrategy,CMAParameters> ipop(func,parameters);
	ipop.optimize();
	return ipop._solutions;
      }
    return CMASolutions();
  };
}

#endif
