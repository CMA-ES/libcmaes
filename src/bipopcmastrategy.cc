
#include "bipopcmastrategy.h"
#include "opti_err.h"
#include <random>
#include <glog/logging.h>
#include <time.h>

namespace libcmaes
{
  template <class TCovarianceUpdate>
  BIPOPCMAStrategy<TCovarianceUpdate>::BIPOPCMAStrategy(FitFunc &func,
							CMAParameters &parameters)
    :IPOPCMAStrategy<TCovarianceUpdate>(func,parameters),_lambda_def(parameters._lambda),_lambda_l(parameters._lambda)
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _gen.seed(static_cast<uint64_t>(time(NULL)));
    _unif = std::uniform_real_distribution<>(0,1);
    //_lambda_def = 4.0+ceil(3.0+log(CMAStrategy<TCovarianceUpdate>::_parameters._dim));
    CMAStrategy<TCovarianceUpdate>::_parameters._lambda = _lambda_def;
    CMAStrategy<TCovarianceUpdate>::_parameters._mu = floor(_lambda_def / 2.0);
    CMAStrategy<TCovarianceUpdate>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate>::_parameters);
  }

  template <class TCovarianceUpdate>
  BIPOPCMAStrategy<TCovarianceUpdate>::~BIPOPCMAStrategy()
  {
  }

  template <class TCovarianceUpdate>
  void BIPOPCMAStrategy<TCovarianceUpdate>::tell()
  {
    CMAStrategy<TCovarianceUpdate>::tell();
  }

  template <class TCovarianceUpdate>
  int BIPOPCMAStrategy<TCovarianceUpdate>::optimize()
  {
    std::array<int,2> budgets = {0,0}; // 0: r1, 1: r2
    CMASolutions best_run;
    for (int r=0;r<CMAStrategy<TCovarianceUpdate>::_parameters._nrestarts;r++)
      {
	while(budgets[0]>budgets[1])
	  {
	    r2();
	    IPOPCMAStrategy<TCovarianceUpdate>::reset_search_state();
	    CMAStrategy<TCovarianceUpdate>::optimize();
	    budgets[1] += CMAStrategy<TCovarianceUpdate>::_solutions._niter * CMAStrategy<TCovarianceUpdate>::_parameters._lambda;
	    IPOPCMAStrategy<TCovarianceUpdate>::capture_best_solution(best_run);
	  }
	if (r > 0) // use lambda_def on first call.
	  {
	    r1();
	    IPOPCMAStrategy<TCovarianceUpdate>::reset_search_state();
	  }
	CMAStrategy<TCovarianceUpdate>::optimize();
	budgets[0] += CMAStrategy<TCovarianceUpdate>::_solutions._niter * CMAStrategy<TCovarianceUpdate>::_parameters._lambda;
	IPOPCMAStrategy<TCovarianceUpdate>::capture_best_solution(best_run);
      }
    CMAStrategy<TCovarianceUpdate>::_solutions = best_run;
    if (CMAStrategy<TCovarianceUpdate>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in CMAStrategy<TCovarianceUpdate>::_solutions._run_status.
  }

  template <class TCovarianceUpdate>
  void BIPOPCMAStrategy<TCovarianceUpdate>::r1()
  {
    CMAStrategy<TCovarianceUpdate>::_parameters._lambda = _lambda_l;
    IPOPCMAStrategy<TCovarianceUpdate>::lambda_inc();
    _lambda_l = CMAStrategy<TCovarianceUpdate>::_parameters._lambda;
  }

  template <class TCovarianceUpdate>
  void BIPOPCMAStrategy<TCovarianceUpdate>::r2()
  {
    double u = _unif(_gen);
    double ltmp = pow(0.5*(_lambda_l/_lambda_def),u);
    double nlambda = ceil(_lambda_def * ltmp);
    LOG_IF(INFO,!CMAStrategy<TCovarianceUpdate>::_parameters._quiet) << "Restart => lambda_s=" << nlambda << " / lambda_old=" << CMAStrategy<TCovarianceUpdate>::_parameters._lambda << " / lambda_l=" << _lambda_l << " / lambda_def=" << _lambda_def << std::endl;
    CMAStrategy<TCovarianceUpdate>::_parameters._lambda = nlambda;
  }

  template class BIPOPCMAStrategy<CovarianceUpdate>;
  template class BIPOPCMAStrategy<ACovarianceUpdate>;
}
