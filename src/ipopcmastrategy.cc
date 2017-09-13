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

#include "ipopcmastrategy.h"
#include "opti_err.h"
#include "llogging.h"
#include <iostream>

namespace libcmaes
{
  template <class TCovarianceUpdate, class TGenoPheno>
  IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::IPOPCMAStrategy(FitFunc &func,
								 CMAParameters<TGenoPheno> &parameters)
    :CMAStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters)
  {
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::IPOPCMAStrategy(FitFunc &func,
								 CMAParameters<TGenoPheno> &parameters,
								 const CMASolutions &solutions)
    :CMAStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters,solutions)
  {
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::~IPOPCMAStrategy()
  {
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::tell()
  {
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::tell();
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  int IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(const EvalFunc &evalf,
							      const AskFunc &askf,
							      const TellFunc &tellf)
  {
    CMASolutions best_run;
    for (int r=0;r<CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._nrestarts;r++)
      {
	LOG_IF(INFO,!(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._quiet)) << "r: " << r << " / lambda=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda << std::endl;
	CMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(evalf,askf,tellf);
		
	// capture best solution.
	capture_best_solution(best_run);
	
	// reset parameters and solutions.
	lambda_inc();
	reset_search_state();
	
	// do not restart if max budget function calls is reached.
	if (CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._max_fevals > 0
	    && CMAStrategy<TCovarianceUpdate,TGenoPheno>::_nevals >= CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._max_fevals)
	  {
	    LOG_IF(INFO,!(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._quiet)) << "IPOP restarts ended on max fevals=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_nevals << ">=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._max_fevals << std::endl;
	    break;
	  }
      }
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = best_run;
    if (CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    return OPTI_ERR_TERMINATION; // exact termination code is in CMAStrategy<TCovarianceUpdate>::_solutions._run_status.
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::lambda_inc()
  {
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda *= 2.0;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters.initialize_parameters();
    LOG_IF(INFO,!(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._quiet)) << "Restart => lambda_l=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda << " / lambda_old=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda / 2.0 << std::endl;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::reset_search_state()
  {
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters);
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_niter = 0;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::capture_best_solution(CMASolutions &best_run)
  {
    if (best_run._candidates.empty() || CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions.best_candidate().get_fvalue() < best_run.best_candidate().get_fvalue())
      best_run = CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions;
  }

  template class CMAES_EXPORT IPOPCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class CMAES_EXPORT IPOPCMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
