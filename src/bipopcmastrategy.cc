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

#include "bipopcmastrategy.h"
#include "opti_err.h"
#include <random>
#include "llogging.h"
#include <ctime>
#include <array>

namespace libcmaes
{
  template <class TCovarianceUpdate, class TGenoPheno>
  BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::BIPOPCMAStrategy(FitFunc &func,
								   CMAParameters<TGenoPheno> &parameters)
    :IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters),_lambda_def(parameters._lambda),_lambda_l(parameters._lambda)
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _gen.seed(static_cast<uint64_t>(time(nullptr)));
    _unif = std::uniform_real_distribution<>(0,1);
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda = _lambda_def;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._mu = floor(_lambda_def / 2.0);
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters);
    _sigma_init = parameters._sigma_init;
    _max_fevals = parameters._max_fevals;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::BIPOPCMAStrategy(FitFunc &func,
								   CMAParameters<TGenoPheno> &parameters,
								   const CMASolutions &solutions)
    :IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters,solutions),_lambda_def(parameters._lambda),_lambda_l(parameters._lambda)
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _gen.seed(static_cast<uint64_t>(time(nullptr)));
    _unif = std::uniform_real_distribution<>(0,1);
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda = _lambda_def;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._mu = floor(_lambda_def / 2.0);
    //CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters);
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::~BIPOPCMAStrategy()
  {
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::tell()
  {
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::tell();
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  int BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(const EvalFunc &evalf,
							       const AskFunc &askf,
							       const TellFunc &tellf)
  {
    std::array<int,2> budgets = {{0,0}}; // 0: r1, 1: r2
    CMASolutions best_run;
    for (int r=0;r<CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._nrestarts;r++)
      {
	while(budgets[0]>budgets[1])
	  {
	    r2();
	    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters.set_max_fevals(0.5*budgets[0]);
	    IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::reset_search_state();
	    CMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(evalf,askf,tellf);
	    budgets[1] += CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions._niter * CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda;
	    IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::capture_best_solution(best_run);
	  }
	if (r > 0) // use lambda_def on first call.
	  {
	    r1();
	    IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::reset_search_state();
	  }
	CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters.set_max_fevals(_max_fevals); // resets the budget
	CMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(evalf,askf,tellf);
	budgets[0] += CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions._niter * CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda;
	IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::capture_best_solution(best_run);
      }
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = best_run;
    if (CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in CMAStrategy<TCovarianceUpdate>::_solutions._run_status.
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::r1()
  {
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda = _lambda_l;
    IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::lambda_inc();
    _lambda_l = CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._sigma_init = _sigma_init;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::r2()
  {
    double u = _unif(_gen);
    double us = _unif(_gen);
    double nsigma = 2.0*pow(10,-2.0*us);
    double ltmp = pow(0.5*(_lambda_l/_lambda_def),u);
    double nlambda = ceil(_lambda_def * ltmp);
    LOG_IF(INFO,!(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._quiet)) << "Restart => lambda_s=" << nlambda << " / lambda_old=" << CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda << " / lambda_l=" << _lambda_l << " / lambda_def=" << _lambda_def << " / nsigma=" << nsigma << std::endl;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._lambda = nlambda;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters._sigma_init = nsigma;
    CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters.initialize_parameters();
  }

  template class CMAES_EXPORT BIPOPCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy> >;
  template class CMAES_EXPORT BIPOPCMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy> >;
}
