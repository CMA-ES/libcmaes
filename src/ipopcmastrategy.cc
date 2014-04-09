
#include "ipopcmastrategy.h"
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{
  IPOPCMAStrategy::IPOPCMAStrategy(FitFunc &func,
				   CMAParameters &parameters)
    :CMAStrategy(func,parameters)
  {
  }

  IPOPCMAStrategy::~IPOPCMAStrategy()
  {
  }
  
  void IPOPCMAStrategy::tell()
  {
    CMAStrategy::tell();
  }

  bool IPOPCMAStrategy::optimize()
  {
    for (int r=0;r<_parameters._nrestarts;r++)
      {
	LOG_IF(INFO,!_parameters._quiet) << "r: " << r << " / lambda=" << _parameters._lambda << std::endl;
	CMAStrategy::optimize();

	// reset parameters and solutions.
	_parameters._lambda *= 2.0;
	_solutions = CMASolutions(_parameters);
	_niter = 0;

	// do not restart if max budget function calls is reached (TODO: or fitness... i.e. if we know the function).
	if (_parameters._max_fevals > 0
	    && _nevals >= _parameters._max_fevals)
	  {
	    LOG_IF(INFO,!_parameters._quiet) << "IPOP restarts ended on max fevals=" << _nevals << ">=" << _parameters._max_fevals << std::endl;
	    break;
	  }
      }
    return true;
  }
}
