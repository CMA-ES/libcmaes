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

#include "bipopcmastrategy.h"
#include "opti_err.h"
#include <random>
#include <glog/logging.h>
#include <time.h>

namespace libcmaes
{
  template <class TCovarianceUpdate, class TBoundStrategy>
  BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::BIPOPCMAStrategy(FitFunc &func,
								       CMAParameters<TBoundStrategy> &parameters)
    :IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>(func,parameters),_lambda_def(parameters._lambda),_lambda_l(parameters._lambda)
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _gen.seed(static_cast<uint64_t>(time(NULL)));
    _unif = std::uniform_real_distribution<>(0,1);
    //_lambda_def = 4.0+ceil(3.0+log(CMAStrategy<TCovarianceUpdate>::_parameters._dim));
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda = _lambda_def;
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._mu = floor(_lambda_def / 2.0);
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate>::_parameters);
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::~BIPOPCMAStrategy()
  {
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  void BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::tell()
  {
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::tell();
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  int BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::optimize()
  {
    std::array<int,2> budgets = {0,0}; // 0: r1, 1: r2
    CMASolutions best_run;
    for (int r=0;r<CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._nrestarts;r++)
      {
	while(budgets[0]>budgets[1])
	  {
	    r2();
	    IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::reset_search_state();
	    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::optimize();
	    budgets[1] += CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_solutions._niter * CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda;
	    IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::capture_best_solution(best_run);
	  }
	if (r > 0) // use lambda_def on first call.
	  {
	    r1();
	    IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::reset_search_state();
	  }
	CMAStrategy<TCovarianceUpdate,TBoundStrategy>::optimize();
	budgets[0] += CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_solutions._niter * CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda;
	IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::capture_best_solution(best_run);
      }
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_solutions = best_run;
    if (CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in CMAStrategy<TCovarianceUpdate>::_solutions._run_status.
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  void BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::r1()
  {
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda = _lambda_l;
    IPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::lambda_inc();
    _lambda_l = CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda;
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  void BIPOPCMAStrategy<TCovarianceUpdate,TBoundStrategy>::r2()
  {
    double u = _unif(_gen);
    double ltmp = pow(0.5*(_lambda_l/_lambda_def),u);
    double nlambda = ceil(_lambda_def * ltmp);
    LOG_IF(INFO,!(CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._quiet)) << "Restart => lambda_s=" << nlambda << " / lambda_old=" << CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda << " / lambda_l=" << _lambda_l << " / lambda_def=" << _lambda_def << std::endl;
    CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_parameters._lambda = nlambda;
  }

  template class BIPOPCMAStrategy<CovarianceUpdate,NoBoundStrategy>;
  template class BIPOPCMAStrategy<ACovarianceUpdate,NoBoundStrategy>;
  //TODO: pwq bound strategy.
}
