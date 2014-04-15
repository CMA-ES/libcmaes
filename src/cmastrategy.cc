
//#define NDEBUG 1

#include "cmastrategy.h"
#include "opti_err.h"
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{
  
  template <class TCovarianceUpdate>
  ProgressFunc<CMAParameters,CMASolutions> CMAStrategy<TCovarianceUpdate>::_defaultPFunc = [](const CMAParameters &cmaparams, const CMASolutions &cmasols)
  {
    LOG_IF(INFO,!cmaparams._quiet) << "iter=" << cmasols._niter << " / evals=" << cmaparams._lambda * cmasols._niter << " / f-value=" << cmasols._best_candidates_hist.back()._fvalue <<  " / sigma=" << cmasols._sigma << (cmaparams._lazy_update && cmasols._updated_eigen ? " / cupdate="+std::to_string(cmasols._updated_eigen) : "");
    return 0;
  };
  
  template <class TCovarianceUpdate>
  CMAStrategy<TCovarianceUpdate>::CMAStrategy(FitFunc &func,
					      CMAParameters &parameters)
    :ESOStrategy(func,parameters)
  {
    _pfunc = _defaultPFunc;
    _esolver = EigenMultivariateNormal<double>(false,_parameters._seed); // seeding the multivariate normal generator.
    LOG_IF(INFO,!_parameters._quiet) << "CMA-ES / dim=" << _parameters._dim << " / lambda=" << _parameters._lambda << " / mu=" << _parameters._mu << " / mueff=" << _parameters._muw << " / c1=" << _parameters._c1 << " / cmu=" << _parameters._cmu << " / lazy_update=" << _parameters._lazy_update << std::endl;
    if (!_parameters._fplot.empty())
      _fplotstream.open(_parameters._fplot);
  }

  template <class TCovarianceUpdate>
  CMAStrategy<TCovarianceUpdate>::~CMAStrategy()
  {
    if (!_parameters._fplot.empty())
      _fplotstream.close();
  }
  
  template <class TCovarianceUpdate>
  dMat CMAStrategy<TCovarianceUpdate>::ask()
  {
    // compute eigenvalues and eigenvectors.
    _solutions._updated_eigen = false;
    if (_niter == 0 || !_parameters._lazy_update
	|| _niter - _solutions._eigeniter > _parameters._lazy_value)
      {
	_solutions._eigeniter = _niter;
	_esolver.setMean(_solutions._xmean);
	_esolver.setCovar(_solutions._cov);
	_solutions._updated_eigen = true;
      }
    
    // sample for multivariate normal distribution.
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
    _solutions.sort_candidates();

    //TODO: test for flat values (same value almost everywhere).

    //TODO: update function value history, as needed.
    _solutions.update_best_candidates();

    //TODO: update best value, as needed.

    // CMA-ES update, depends on the selected 'flavor'.
    TCovarianceUpdate::update(_parameters,_esolver,_solutions);
    
    // other stuff.
    _solutions.update_eigenv(_esolver._eigenSolver.eigenvalues(),
			     _esolver._eigenSolver.eigenvectors());
    _solutions._niter = _niter;
  }

  template <class TCovarianceUpdate>
  bool CMAStrategy<TCovarianceUpdate>::stop()
  {
    if (_solutions._run_status < 0) // an error occured, most likely out of memory at cov matrix creation.
      return true;
    
    if (_niter == 0)
      return false;
    
    if (_pfunc(_parameters,_solutions)) // progress function.
      return true; // end on progress function internal termination, possibly custom.
    
    if (!_parameters._fplot.empty())
      plot();
    
    if ((_parameters._max_iter > 0 && _niter >= _parameters._max_iter)
	|| (_solutions._run_status = _stopcriteria.stop(_parameters,_solutions)) != 0)
      return true;
    else return false;
  }

  template <class TCovarianceUpdate>
  int CMAStrategy<TCovarianceUpdate>::optimize()
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
    if (_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in _solutions._run_status.
  }

  template <class TCovarianceUpdate>
  void CMAStrategy<TCovarianceUpdate>::plot()
  {
    static std::string sep = " ";
    _fplotstream << fabs(_solutions._best_candidates_hist.back()._fvalue) << sep
		 << _nevals << sep << _solutions._sigma << sep << sqrt(_solutions._max_eigenv/_solutions._min_eigenv) << sep;
    _fplotstream << _esolver._eigenSolver.eigenvalues().transpose() << sep; // eigenvalues
    _fplotstream << _solutions._sigma * _solutions._cov.colwise().maxCoeff().array().sqrt() << sep; // max deviation in all main axes
    _fplotstream << _solutions._xmean.transpose();
    _fplotstream << std::endl;
  }
  
  template class CMAStrategy<CovarianceUpdate>;
  
}
