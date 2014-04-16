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
