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

//#define NDEBUG 1

#include "cmastrategy.h"
#include "opti_err.h"
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{

  template <class TBoundStrategy> using eostrat = ESOStrategy<CMAParameters<TBoundStrategy>,CMASolutions,CMAStopCriteria<TBoundStrategy> >;
  
  template <class TCovarianceUpdate, class TBoundStrategy>
  ProgressFunc<CMAParameters<TBoundStrategy>,CMASolutions> CMAStrategy<TCovarianceUpdate,TBoundStrategy>::_defaultPFunc = [](const CMAParameters<TBoundStrategy> &cmaparams, const CMASolutions &cmasols)
  {
    LOG_IF(INFO,!cmaparams._quiet) << "iter=" << cmasols._niter << " / evals=" << cmaparams._lambda * cmasols._niter << " / f-value=" << cmasols._best_candidates_hist.back()._fvalue <<  " / sigma=" << cmasols._sigma << (cmaparams._lazy_update && cmasols._updated_eigen ? " / cupdate="+std::to_string(cmasols._updated_eigen) : "");
    return 0;
  };
  
  template <class TCovarianceUpdate, class TBoundStrategy>
  CMAStrategy<TCovarianceUpdate,TBoundStrategy>::CMAStrategy(FitFunc &func,
					      CMAParameters<TBoundStrategy> &parameters)
    :ESOStrategy<CMAParameters<TBoundStrategy>,CMASolutions,CMAStopCriteria<TBoundStrategy> >(func,parameters)
  {
    eostrat<TBoundStrategy>::_pfunc = _defaultPFunc;
    _esolver = EigenMultivariateNormal<double>(false,eostrat<TBoundStrategy>::_parameters._seed); // seeding the multivariate normal generator.
    LOG_IF(INFO,!eostrat<TBoundStrategy>::_parameters._quiet) << "CMA-ES / dim=" << eostrat<TBoundStrategy>::_parameters._dim << " / lambda=" << eostrat<TBoundStrategy>::_parameters._lambda << " / mu=" << eostrat<TBoundStrategy>::_parameters._mu << " / mueff=" << eostrat<TBoundStrategy>::_parameters._muw << " / c1=" << eostrat<TBoundStrategy>::_parameters._c1 << " / cmu=" << eostrat<TBoundStrategy>::_parameters._cmu << " / lazy_update=" << eostrat<TBoundStrategy>::_parameters._lazy_update << std::endl;
    if (!eostrat<TBoundStrategy>::_parameters._fplot.empty())
      _fplotstream.open(eostrat<TBoundStrategy>::_parameters._fplot);
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  CMAStrategy<TCovarianceUpdate,TBoundStrategy>::~CMAStrategy()
  {
    if (!eostrat<TBoundStrategy>::_parameters._fplot.empty())
      _fplotstream.close();
  }
  
  template <class TCovarianceUpdate, class TBoundStrategy>
  dMat CMAStrategy<TCovarianceUpdate,TBoundStrategy>::ask()
  {
    // compute eigenvalues and eigenvectors.
    eostrat<TBoundStrategy>::_solutions._updated_eigen = false;
    if (eostrat<TBoundStrategy>::_niter == 0 || !eostrat<TBoundStrategy>::_parameters._lazy_update
	|| eostrat<TBoundStrategy>::_niter - eostrat<TBoundStrategy>::_solutions._eigeniter > eostrat<TBoundStrategy>::_parameters._lazy_value)
      {
	eostrat<TBoundStrategy>::_solutions._eigeniter = eostrat<TBoundStrategy>::_niter;
	_esolver.setMean(eostrat<TBoundStrategy>::_solutions._xmean);
	_esolver.setCovar(eostrat<TBoundStrategy>::_solutions._cov);
	eostrat<TBoundStrategy>::_solutions._updated_eigen = true;
      }
    
    // sample for multivariate normal distribution.
    dMat pop = _esolver.samples(eostrat<TBoundStrategy>::_parameters._lambda,eostrat<TBoundStrategy>::_solutions._sigma); // Eq (1).
    
    //TODO: rescale to function space as needed.

    //debug
    /*DLOG(INFO) << "ask: produced " << pop.cols() << " candidates\n";
      std::cerr << pop << std::endl;*/
    //debug
    
    return pop;
  }
  
  template <class TCovarianceUpdate, class TBoundStrategy>
  void CMAStrategy<TCovarianceUpdate,TBoundStrategy>::tell()
  {
    //debug
    //DLOG(INFO) << "tell()\n";
    //debug
    
    // sort candidates.
    eostrat<TBoundStrategy>::_solutions.sort_candidates();

    //TODO: test for flat values (same value almost everywhere).

    //TODO: update function value history, as needed.
    eostrat<TBoundStrategy>::_solutions.update_best_candidates();

    //TODO: update best value, as needed.

    // CMA-ES update, depends on the selected 'flavor'.
    TCovarianceUpdate::update(eostrat<TBoundStrategy>::_parameters,_esolver,eostrat<TBoundStrategy>::_solutions);
    
    // other stuff.
    eostrat<TBoundStrategy>::_solutions.update_eigenv(_esolver._eigenSolver.eigenvalues(),
			     _esolver._eigenSolver.eigenvectors());
    eostrat<TBoundStrategy>::_solutions._niter = eostrat<TBoundStrategy>::_niter;
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  bool CMAStrategy<TCovarianceUpdate,TBoundStrategy>::stop()
  {
    if (eostrat<TBoundStrategy>::_solutions._run_status < 0) // an error occured, most likely out of memory at cov matrix creation.
      return true;
    
    if (eostrat<TBoundStrategy>::_niter == 0)
      return false;
    
    if (eostrat<TBoundStrategy>::_pfunc(eostrat<TBoundStrategy>::_parameters,eostrat<TBoundStrategy>::_solutions)) // progress function.
      return true; // end on progress function internal termination, possibly custom.
    
    if (!eostrat<TBoundStrategy>::_parameters._fplot.empty())
      plot();
    
    if ((eostrat<TBoundStrategy>::_parameters._max_iter > 0 && eostrat<TBoundStrategy>::_niter >= eostrat<TBoundStrategy>::_parameters._max_iter)
	|| (eostrat<TBoundStrategy>::_solutions._run_status = _stopcriteria.stop(eostrat<TBoundStrategy>::_parameters,eostrat<TBoundStrategy>::_solutions)) != 0)
      return true;
    else return false;
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  int CMAStrategy<TCovarianceUpdate,TBoundStrategy>::optimize()
  {
    //debug
    //DLOG(INFO) << "optimize()\n";
    //debug
    
    while(!stop())
      {
	dMat candidates = ask();
	this->eval(eostrat<TBoundStrategy>::_parameters._gp.pheno(candidates));
	tell();
	eostrat<TBoundStrategy>::_niter++;
      }
    if (eostrat<TBoundStrategy>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in eostrat<TBoundStrategy>::_solutions._run_status.
  }

  template <class TCovarianceUpdate, class TBoundStrategy>
  void CMAStrategy<TCovarianceUpdate,TBoundStrategy>::plot()
  {
    static std::string sep = " ";
    _fplotstream << fabs(eostrat<TBoundStrategy>::_solutions._best_candidates_hist.back()._fvalue) << sep
		 << eostrat<TBoundStrategy>::_nevals << sep << eostrat<TBoundStrategy>::_solutions._sigma << sep << sqrt(eostrat<TBoundStrategy>::_solutions._max_eigenv/eostrat<TBoundStrategy>::_solutions._min_eigenv) << sep;
    _fplotstream << _esolver._eigenSolver.eigenvalues().transpose() << sep; // eigenvalues
    _fplotstream << eostrat<TBoundStrategy>::_solutions._cov.colwise().maxCoeff().array().sqrt() << sep; // max deviation in all main axes
    _fplotstream << eostrat<TBoundStrategy>::_solutions._xmean.transpose();
    _fplotstream << std::endl;
  }
  
  template class CMAStrategy<CovarianceUpdate,NoBoundStrategy>;
  template class CMAStrategy<ACovarianceUpdate,NoBoundStrategy>;
  //TODO: pwq bound strategy.
}
