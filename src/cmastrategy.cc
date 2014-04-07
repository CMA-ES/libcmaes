
//#define NDEBUG 1

#include "cmastrategy.h"
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{

  template <class TCovarianceUpdate>
  CMAStrategy<TCovarianceUpdate>::CMAStrategy(FitFunc &func,
					      CMAParameters &parameters)
    :ESOStrategy(func,parameters)
  {    
    LOG_IF(INFO,!_parameters._quiet) << "CMA-ES / dim=" << _parameters._dim << " / lambda=" << _parameters._lambda << " / mu=" << _parameters._mu << " / mueff=" << _parameters._muw << std::endl;
  }

  template <class TCovarianceUpdate>
  CMAStrategy<TCovarianceUpdate>::~CMAStrategy()
  {
  }
  
  template <class TCovarianceUpdate>
  dMat CMAStrategy<TCovarianceUpdate>::ask()
  {
    // sample for multivariate normal distribution.
    _esolver = EigenMultivariateNormal<double>(_solutions._xmean,_solutions._cov);
    dMat pop = _esolver.samples(_parameters._lambda,_solutions._sigma); // Eq (1).
    
    //TODO: rescale to function space as needed.

    //debug
    /*DLOG(INFO) << "ask: produced " << pop.cols() << " candidates\n";
      std::cerr << pop << std::endl;*/
    //debug
    
    return pop;
  }
  
  template <class TCovarianceUpdate>
  void CMAStrategy<TCovarianceUpdate>::tell()
  {
    //debug
    //DLOG(INFO) << "tell()\n";
    //debug
    
    // sort candidates.
    _solutions.sort();

    //TODO: test for flat values (same value almost everywhere).

    //TODO: update function value history, as needed.
    _solutions.update_best_candidates();

    //TODO: update best value, as needed.

    // compute mean, Eq. (2)
    dVec xmean = dVec::Zero(_parameters._dim);
    for (int i=0;i<_parameters._mu;i++)
      xmean += _parameters._weights[i] * _solutions._candidates.at(i)._x;
    
    // reusable variables.
    dVec diffxmean = 1.0/_solutions._sigma * (xmean-_solutions._xmean); // (m^{t+1}-m^t)/sigma^t
    dMat Csqinv = _esolver._eigenSolver.operatorInverseSqrt();

    // update psigma, Eq. (3)
    _solutions._psigma = (1.0-_parameters._csigma)*_solutions._psigma
      + sqrt(_parameters._csigma*(2.0-_parameters._csigma)*_parameters._muw) * Csqinv * diffxmean; //TODO: use _facs_ps instead of sqrt...
    double norm_ps = _solutions._psigma.norm();

    // update pc, Eq. (4)
    _solutions._hsig = 1;
    double val_for_hsig = sqrt(1.0-pow(1.0-_parameters._csigma,2.0*(_niter+1)))*(1.4+2.0/(_parameters._dim+1))*_parameters._chi;
    if (norm_ps < val_for_hsig)
      _solutions._hsig = 0; //TODO: simplify equation instead.
    _solutions._pc = (1-_parameters._cc) * _solutions._pc + _solutions._hsig * sqrt(_parameters._cc * (2.0 - _parameters._cc) * _parameters._muw) * diffxmean;
    dMat spc = _solutions._pc * _solutions._pc.transpose();
    
    // covariance update, Eq (5).
    dMat wdiff = dMat::Zero(_parameters._dim,_parameters._dim);
    for (int i=0;i<_parameters._mu;i++)
      {
	dVec difftmp = _solutions._candidates.at(i)._x - _solutions._xmean;
	wdiff += _parameters._weights[i] * (difftmp*difftmp.transpose());
      }
    wdiff *= 1.0/(_solutions._sigma*_solutions._sigma);
    _solutions._cov = (1-_parameters._c1-_parameters._cmu+(1-_solutions._hsig)*_parameters._c1*_parameters._cc*(2.0-_parameters._cc))*_solutions._cov + _parameters._c1*spc + _parameters._cmu*wdiff;
    
    // sigma update, Eq. (6)
    _solutions._sigma *= exp((_parameters._csigma / _parameters._dsigma) * (norm_ps / _parameters._chi - 1.0));
    
    // set mean.
    _solutions._xmean = xmean;
  }

  template <class TCovarianceUpdate>
  bool CMAStrategy<TCovarianceUpdate>::stop()
  {
    if (_niter == 0)
      return false;
    
    //TODO: other termination criterions.
    double diff_value = 100.0;
    if (_solutions._best_candidates_hist.size()>1)
      {
	diff_value = _solutions._best_candidates_hist.back()._fvalue - _solutions._best_candidates_hist.at(_solutions._best_candidates_hist.size()-2)._fvalue;
      }

    LOG_IF(INFO,!_parameters._quiet) << "iter=" << _niter << " / evals=" << _nevals << " / f-value=" << _solutions._best_candidates_hist.back()._fvalue <<  " / sigma=" << _solutions._sigma << " / f-diff value=" << diff_value << std::endl;
    
    // for now a simple diff in func value, for testing purposes.
    if (/*fabs(diff_value) < 1e-12
	  ||*/ (_parameters._max_iter > 0 && _niter >= _parameters._max_iter))
      return true;
    else return false;
  }

  template <class TCovarianceUpdate>
  bool CMAStrategy<TCovarianceUpdate>::optimize()
  {
    //debug
    //DLOG(INFO) << "optimize()\n";
    //debug
    
    while(!stop())
      {
	dMat candidates = ask();
	eval(candidates);
	tell();
	_niter++;
      }
    return true;
  }
  
  template class CMAStrategy<CovarianceUpdate>;
  
}
