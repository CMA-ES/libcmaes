
#include "ipopcmastrategy.h"
#include "opti_err.h"
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{
  template <class TCovarianceUpdate>
  IPOPCMAStrategy<TCovarianceUpdate>::IPOPCMAStrategy(FitFunc &func,
						      CMAParameters &parameters)
    :CMAStrategy<TCovarianceUpdate>(func,parameters)
  {
  }

  template <class TCovarianceUpdate>
  IPOPCMAStrategy<TCovarianceUpdate>::~IPOPCMAStrategy()
  {
  }

  template <class TCovarianceUpdate>
  void IPOPCMAStrategy<TCovarianceUpdate>::tell()
  {
    CMAStrategy<TCovarianceUpdate>::tell();
  }

  template <class TCovarianceUpdate>
  int IPOPCMAStrategy<TCovarianceUpdate>::optimize()
  {
    CMASolutions best_run;
    for (int r=0;r<CMAStrategy<TCovarianceUpdate>::_parameters._nrestarts;r++)
      {
	LOG_IF(INFO,!CMAStrategy<TCovarianceUpdate>::_parameters._quiet) << "r: " << r << " / lambda=" << CMAStrategy<TCovarianceUpdate>::_parameters._lambda << std::endl;
	CMAStrategy<TCovarianceUpdate>::optimize();

	// capture best solution.
	capture_best_solution(best_run);
	
	// reset parameters and solutions.
	lambda_inc();
	reset_search_state();
	
	// do not restart if max budget function calls is reached (TODO: or fitness... i.e. if we know the function).
	if (CMAStrategy<TCovarianceUpdate>::_parameters._max_fevals > 0
	    && CMAStrategy<TCovarianceUpdate>::_nevals >= CMAStrategy<TCovarianceUpdate>::_parameters._max_fevals)
	  {
	    LOG_IF(INFO,!CMAStrategy<TCovarianceUpdate>::_parameters._quiet) << "IPOP restarts ended on max fevals=" << CMAStrategy<TCovarianceUpdate>::_nevals << ">=" << CMAStrategy<TCovarianceUpdate>::_parameters._max_fevals << std::endl;
	    break;
	  }
      }
    CMAStrategy<TCovarianceUpdate>::_solutions = best_run;
    if (CMAStrategy<TCovarianceUpdate>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    return OPTI_ERR_TERMINATION; // exact termination code is in CMAStrategy<TCovarianceUpdate>::_solutions._run_status.
  }

  template <class TCovarianceUpdate>
  void IPOPCMAStrategy<TCovarianceUpdate>::lambda_inc()
  {
    CMAStrategy<TCovarianceUpdate>::_parameters._lambda *= 2.0;
    LOG_IF(INFO,!CMAStrategy<TCovarianceUpdate>::_parameters._quiet) << "Restart => lambda_l=" << CMAStrategy<TCovarianceUpdate>::_parameters._lambda << " / lambda_old=" << CMAStrategy<TCovarianceUpdate>::_parameters._lambda / 2.0 << std::endl;
  }

  template <class TCovarianceUpdate>
  void IPOPCMAStrategy<TCovarianceUpdate>::reset_search_state()
  {
    CMAStrategy<TCovarianceUpdate>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate>::_parameters);
    CMAStrategy<TCovarianceUpdate>::_niter = 0;
  }

  template <class TCovarianceUpdate>
  void IPOPCMAStrategy<TCovarianceUpdate>::capture_best_solution(CMASolutions &best_run)
  {
    if (best_run._candidates.empty() || CMAStrategy<TCovarianceUpdate>::_solutions.best_candidate()._fvalue < best_run.best_candidate()._fvalue)
      best_run = CMAStrategy<TCovarianceUpdate>::_solutions;
  }

  template class IPOPCMAStrategy<CovarianceUpdate>;
  template class IPOPCMAStrategy<ACovarianceUpdate>;
}
