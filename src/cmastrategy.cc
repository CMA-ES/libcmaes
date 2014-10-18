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

//#define NDEBUG 1

#include "libcmaes_config.h"
#include "cmastrategy.h"
#include "opti_err.h"
#include "llogging.h"
#include <iostream>
#include <chrono>

namespace libcmaes
{

  template <class TGenoPheno> using eostrat = ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >;
  
  template <class TCovarianceUpdate, class TGenoPheno>
  ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> CMAStrategy<TCovarianceUpdate,TGenoPheno>::_defaultPFunc = [](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols)
  {
    LOG_IF(INFO,!cmaparams.quiet()) << "iter=" << cmasols.niter() << " / evals=" << cmasols.fevals() << " / f-value=" << cmasols.best_candidate().get_fvalue() <<  " / sigma=" << cmasols.sigma() << " / last_iter=" << cmasols.elapsed_last_iter() << std::endl;
    return 0;
  };

  template<class TCovarianceUpdate, class TGenoPheno>
  PlotFunc<CMAParameters<TGenoPheno>,CMASolutions> CMAStrategy<TCovarianceUpdate,TGenoPheno>::_defaultFPFunc = [](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
  {
    std::string sep = " ";
    fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << sqrt(cmasols.max_eigenv()/cmasols.min_eigenv()) << sep;
    fplotstream << cmasols.eigenvalues().transpose() << sep;
    if (!cmaparams.is_sep())
      fplotstream << cmasols.cov().sqrt().diagonal().transpose() << sep; // max deviation in all main axes
    else fplotstream << cmasols.sepcov().cwiseSqrt().transpose() << sep;
    fplotstream << cmaparams.get_gp().pheno(cmasols.xmean()).transpose();
    fplotstream << sep << cmasols.elapsed_last_iter();
#ifdef HAVE_DEBUG
    fplotstream << sep << cmasols._elapsed_eval << sep << cmasols._elapsed_ask << sep << cmasols._elapsed_tell << sep << cmasols._elapsed_stop;
#endif
    fplotstream << std::endl;
    return 0;
  };

  template <class TCovarianceUpdate, class TGenoPheno>
  CMAStrategy<TCovarianceUpdate,TGenoPheno>::CMAStrategy()
    :ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >()
  {
  }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  CMAStrategy<TCovarianceUpdate,TGenoPheno>::CMAStrategy(FitFunc &func,
							 CMAParameters<TGenoPheno> &parameters)
    :ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >(func,parameters)
  {
    eostrat<TGenoPheno>::_pfunc = [](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols)
      {
	LOG_IF(INFO,!cmaparams.quiet()) << "iter=" << cmasols.niter() << " / evals=" << cmasols.fevals() << " / f-value=" << cmasols.best_candidate().get_fvalue() <<  " / sigma=" << cmasols.sigma() << " / last_iter=" << cmasols.elapsed_last_iter() << std::endl;
	return 0;
      };
    eostrat<TGenoPheno>::_pffunc = [](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
      {
	std::string sep = " ";
	fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << sqrt(cmasols.max_eigenv()/cmasols.min_eigenv()) << sep;
	fplotstream << cmasols.eigenvalues().transpose() << sep;
	if (!cmaparams.is_sep())
	  fplotstream << cmasols.cov().sqrt().diagonal().transpose() << sep; // max deviation in all main axes
	else fplotstream << cmasols.sepcov().cwiseSqrt().transpose() << sep;
	fplotstream << cmaparams.get_gp().pheno(cmasols.xmean()).transpose();
	fplotstream << sep << cmasols.elapsed_last_iter();
#ifdef HAVE_DEBUG
	fplotstream << sep << cmasols._elapsed_eval << sep << cmasols._elapsed_ask << sep << cmasols._elapsed_tell << sep << cmasols._elapsed_stop;
#endif
	fplotstream << std::endl;
	return 0;
      };
    _esolver = EigenMultivariateNormal<double>(false,eostrat<TGenoPheno>::_parameters._seed); // seeding the multivariate normal generator.
    LOG_IF(INFO,!eostrat<TGenoPheno>::_parameters._quiet) << "CMA-ES / dim=" << eostrat<TGenoPheno>::_parameters._dim << " / lambda=" << eostrat<TGenoPheno>::_parameters._lambda << " / sigma0=" << eostrat<TGenoPheno>::_solutions._sigma << " / mu=" << eostrat<TGenoPheno>::_parameters._mu << " / mueff=" << eostrat<TGenoPheno>::_parameters._muw << " / c1=" << eostrat<TGenoPheno>::_parameters._c1 << " / cmu=" << eostrat<TGenoPheno>::_parameters._cmu << " / lazy_update=" << eostrat<TGenoPheno>::_parameters._lazy_update << std::endl;
    if (!eostrat<TGenoPheno>::_parameters._fplot.empty())
      _fplotstream = new std::ofstream(eostrat<TGenoPheno>::_parameters._fplot);
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  CMAStrategy<TCovarianceUpdate,TGenoPheno>::~CMAStrategy()
  {
    if (!eostrat<TGenoPheno>::_parameters._fplot.empty())
      delete _fplotstream;
  }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  dMat CMAStrategy<TCovarianceUpdate,TGenoPheno>::ask()
  {
#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
#endif
    
    // compute eigenvalues and eigenvectors.
    if (!eostrat<TGenoPheno>::_parameters._sep)
      {
	eostrat<TGenoPheno>::_solutions._updated_eigen = false;
	if (eostrat<TGenoPheno>::_niter == 0 || !eostrat<TGenoPheno>::_parameters._lazy_update
	    || eostrat<TGenoPheno>::_niter - eostrat<TGenoPheno>::_solutions._eigeniter > eostrat<TGenoPheno>::_parameters._lazy_value)
	  {
	    eostrat<TGenoPheno>::_solutions._eigeniter = eostrat<TGenoPheno>::_niter;
	    _esolver.setMean(eostrat<TGenoPheno>::_solutions._xmean);
	    _esolver.setCovar(eostrat<TGenoPheno>::_solutions._cov);
	    eostrat<TGenoPheno>::_solutions._updated_eigen = true;
	  }
      }
    else
      {
	_esolver.setMean(eostrat<TGenoPheno>::_solutions._xmean);
	_esolver.set_covar(eostrat<TGenoPheno>::_solutions._sepcov);
	_esolver.set_transform(eostrat<TGenoPheno>::_solutions._sepcov.cwiseSqrt());
      }

    //debug
    //std::cout << "transform: " << _esolver._transform << std::endl;
    //debug
    
    // sample for multivariate normal distribution, produces one candidate per column.
    dMat pop;
    if (!eostrat<TGenoPheno>::_parameters._sep)
      pop = _esolver.samples(eostrat<TGenoPheno>::_parameters._lambda,eostrat<TGenoPheno>::_solutions._sigma); // Eq (1).
    else pop = _esolver.samples_ind(eostrat<TGenoPheno>::_parameters._lambda,eostrat<TGenoPheno>::_solutions._sigma);

    // gradient if available.
    if (eostrat<TGenoPheno>::_parameters._with_gradient)
      {
	dVec grad_at_mean = eostrat<TGenoPheno>::gradf(eostrat<TGenoPheno>::_parameters._gp.pheno(eostrat<TGenoPheno>::_solutions._xmean));
	dVec gradgp_at_mean = eostrat<TGenoPheno>::gradgp(eostrat<TGenoPheno>::_solutions._xmean); // for geno / pheno transform.
	grad_at_mean = grad_at_mean.cwiseProduct(gradgp_at_mean);
	if (grad_at_mean != dVec::Zero(eostrat<TGenoPheno>::_parameters._dim))
	  {
	    dVec nx;
	    if (!eostrat<TGenoPheno>::_parameters._sep)
	      {
		dMat sqrtcov = _esolver._eigenSolver.operatorSqrt();
		dVec q = sqrtcov * grad_at_mean;
		double normq = q.squaredNorm();
		nx = eostrat<TGenoPheno>::_solutions._xmean - eostrat<TGenoPheno>::_solutions._sigma * (sqrt(eostrat<TGenoPheno>::_parameters._dim / normq)) * eostrat<TGenoPheno>::_solutions._cov * grad_at_mean;
	      }
	    else nx = eostrat<TGenoPheno>::_solutions._xmean - eostrat<TGenoPheno>::_solutions._sigma * (sqrt(eostrat<TGenoPheno>::_parameters._dim) / ((eostrat<TGenoPheno>::_solutions._sepcov.cwiseSqrt().cwiseProduct(grad_at_mean)).norm())) * eostrat<TGenoPheno>::_solutions._sepcov.cwiseProduct(grad_at_mean);
	    pop.col(0) = nx;
	  }
      }
    
    // if some parameters are fixed, reset them.
    if (!eostrat<TGenoPheno>::_parameters._fixed_p.empty())
      {
	for (auto it=eostrat<TGenoPheno>::_parameters._fixed_p.begin();
	     it!=eostrat<TGenoPheno>::_parameters._fixed_p.end();++it)
	  {
	    pop.block((*it).first,0,1,pop.cols()) = dVec::Constant(pop.cols(),(*it).second).transpose();
	  }
      }
    
    //debug
    /*DLOG(INFO) << "ask: produced " << pop.cols() << " candidates\n";
      std::cerr << pop << std::endl;*/
    //debug

#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
    eostrat<TGenoPheno>::_solutions._elapsed_ask = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
    
    return pop;
  }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  void CMAStrategy<TCovarianceUpdate,TGenoPheno>::tell()
  {
    //debug
    //DLOG(INFO) << "tell()\n";
    //debug

#ifdef DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
#endif
    
    // sort candidates.
    eostrat<TGenoPheno>::_solutions.sort_candidates();

    // update function value history, as needed.
    eostrat<TGenoPheno>::_solutions.update_best_candidates();
    
    // CMA-ES update, depends on the selected 'flavor'.
    TCovarianceUpdate::update(eostrat<TGenoPheno>::_parameters,_esolver,eostrat<TGenoPheno>::_solutions);
    
    // other stuff.
    if (!eostrat<TGenoPheno>::_parameters._sep)
      eostrat<TGenoPheno>::_solutions.update_eigenv(_esolver._eigenSolver.eigenvalues(),
						    _esolver._eigenSolver.eigenvectors());
    else eostrat<TGenoPheno>::_solutions.update_eigenv(eostrat<TGenoPheno>::_solutions._sepcov,
						       dMat::Constant(eostrat<TGenoPheno>::_parameters._dim,1,1.0));
#ifdef DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
    eostrat<TGenoPheno>::_solutions._elapsed_tell = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  bool CMAStrategy<TCovarianceUpdate,TGenoPheno>::stop()
  {
    if (eostrat<TGenoPheno>::_solutions._run_status < 0) // an error occured, most likely out of memory at cov matrix creation.
      return true;
    
    if (eostrat<TGenoPheno>::_niter == 0)
      return false;
    
    if (eostrat<TGenoPheno>::_pfunc(eostrat<TGenoPheno>::_parameters,eostrat<TGenoPheno>::_solutions)) // progress function.
      return true; // end on progress function internal termination, possibly custom.
    
    if (!eostrat<TGenoPheno>::_parameters._fplot.empty())
      plot();
    
    if ((eostrat<TGenoPheno>::_solutions._run_status = _stopcriteria.stop(eostrat<TGenoPheno>::_parameters,eostrat<TGenoPheno>::_solutions)) != CONT)
      return true;
    else return false;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  int CMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize()
  {
    //debug
    //DLOG(INFO) << "optimize()\n";
    //debug
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
    while(!stop())
      {
	dMat candidates = ask();
	this->eval(candidates,eostrat<TGenoPheno>::_parameters._gp.pheno(candidates));
	tell();
	eostrat<TGenoPheno>::inc_iter();
	std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	eostrat<TGenoPheno>::_solutions._elapsed_last_iter = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	tstart = std::chrono::system_clock::now();
      }
    if (eostrat<TGenoPheno>::_parameters._with_edm)
      eostrat<TGenoPheno>::edm();
    if (eostrat<TGenoPheno>::_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in eostrat<TGenoPheno>::_solutions._run_status.
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  void CMAStrategy<TCovarianceUpdate,TGenoPheno>::plot()
  {
    eostrat<TGenoPheno>::_pffunc(eostrat<TGenoPheno>::_parameters,eostrat<TGenoPheno>::_solutions,*_fplotstream);
  }
  
  template class CMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
