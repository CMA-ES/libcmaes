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


#include "libcmaes_config.h"
#include "cmastrategy.h"
#include "opti_err.h"
#include "llogging.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace libcmaes
{

  template <class TGenoPheno> using eostrat = ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >;
  
  template <class TCovarianceUpdate, class TGenoPheno>
  int pfuncdef_impl(const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols)
  {
    LOG_IF(INFO,!cmaparams.quiet()) << std::setprecision(std::numeric_limits<double>::digits10) << "iter=" << cmasols.niter() << " / evals=" << cmasols.fevals() << " / f-value=" << cmasols.best_candidate().get_fvalue() <<  " / sigma=" << cmasols.sigma() << " / last_iter=" << cmasols.elapsed_last_iter() << std::endl;
    return 0;
  }
  template <class TCovarianceUpdate, class TGenoPheno>
  ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> CMAStrategy<TCovarianceUpdate,TGenoPheno>::_defaultPFunc = &pfuncdef_impl<TCovarianceUpdate,TGenoPheno>;
  
  template<class TCovarianceUpdate, class TGenoPheno>
  int fpfuncdef_full_impl(const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
  {
    std::string sep = " ";
    if (cmasols.niter() == 0)
      {
	std::chrono::time_point<std::chrono::system_clock> tnow = std::chrono::system_clock::now();
	std::time_t tdate = std::chrono::system_clock::to_time_t(tnow);
	fplotstream << cmaparams.dim() << sep << cmaparams.get_seed() << " / " << std::ctime(&tdate) << std::endl; // date and seed
      }
    fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << (cmasols.min_eigenv() == 0 ? 1.0 : sqrt(cmasols.max_eigenv()/cmasols.min_eigenv())) << sep;
    fplotstream << cmasols.get_best_seen_candidate().get_fvalue() << sep << cmasols.get_candidate(cmasols.size() / 2).get_fvalue() << sep << cmasols.get_worst_seen_candidate().get_fvalue() << sep << cmasols.min_eigenv() << sep << cmasols.max_eigenv() << sep; // best ever fvalue, median fvalue, worst fvalue, max eigen, min eigen
    if (cmasols.get_best_seen_candidate().get_x_size())
      fplotstream << cmasols.get_best_seen_candidate().get_x_dvec().transpose() << sep;  // xbest
    else fplotstream << dVec::Zero(cmaparams.dim()).transpose() << sep;
    if (!cmasols.eigenvalues().size())
      fplotstream << dVec::Zero(cmaparams.dim()).transpose() << sep;
    else fplotstream << cmasols.eigenvalues().transpose() << sep;
    fplotstream << cmasols.stds(cmaparams).transpose() << sep;
    fplotstream << cmaparams.get_gp().pheno(cmasols.xmean()).transpose();
    fplotstream << sep << cmasols.elapsed_last_iter();
#ifdef HAVE_DEBUG
    fplotstream << sep << cmasols._elapsed_eval << sep << cmasols._elapsed_ask << sep << cmasols._elapsed_tell << sep << cmasols._elapsed_stop;
#endif
    fplotstream << std::endl;
    return 0;
    }
  template<class TCovarianceUpdate, class TGenoPheno>
  int fpfuncdef_impl(const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
  {
    std::string sep = " ";
    fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << (cmasols.min_eigenv() == 0 ? 1.0 : sqrt(cmasols.max_eigenv()/cmasols.min_eigenv())) << sep;
    if (!cmasols.eigenvalues().size())
      fplotstream << dVec::Zero(cmaparams.dim()).transpose() << sep;
    else fplotstream << cmasols.eigenvalues().transpose() << sep;
    fplotstream << cmasols.stds(cmaparams).transpose() << sep;
    fplotstream << cmaparams.get_gp().pheno(cmasols.xmean()).transpose();
    fplotstream << sep << cmasols.elapsed_last_iter();
#ifdef HAVE_DEBUG
    fplotstream << sep << cmasols._elapsed_eval << sep << cmasols._elapsed_ask << sep << cmasols._elapsed_tell << sep << cmasols._elapsed_stop;
#endif
    fplotstream << std::endl;
    return 0;
    }
  template<class TCovarianceUpdate, class TGenoPheno>
  PlotFunc<CMAParameters<TGenoPheno>,CMASolutions> CMAStrategy<TCovarianceUpdate,TGenoPheno>::_defaultFPFunc = &fpfuncdef_impl<TCovarianceUpdate,TGenoPheno>;

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
    eostrat<TGenoPheno>::_pfunc = _defaultPFunc;
    if (!parameters._full_fplot)
      eostrat<TGenoPheno>::_pffunc = _defaultFPFunc;
    else eostrat<TGenoPheno>::_pffunc = &fpfuncdef_full_impl<TCovarianceUpdate,TGenoPheno>;
    _esolver = Eigen::EigenMultivariateNormal<double>(false,eostrat<TGenoPheno>::_parameters._seed); // seeding the multivariate normal generator.
    LOG_IF(INFO,!eostrat<TGenoPheno>::_parameters._quiet) << "CMA-ES / dim=" << eostrat<TGenoPheno>::_parameters._dim << " / lambda=" << eostrat<TGenoPheno>::_parameters._lambda << " / sigma0=" << eostrat<TGenoPheno>::_solutions._sigma << " / mu=" << eostrat<TGenoPheno>::_parameters._mu << " / mueff=" << eostrat<TGenoPheno>::_parameters._muw << " / c1=" << eostrat<TGenoPheno>::_parameters._c1 << " / cmu=" << eostrat<TGenoPheno>::_parameters._cmu << " / tpa=" << (eostrat<TGenoPheno>::_parameters._tpa==2) << " / threads=" << Eigen::nbThreads() << std::endl;
    if (!eostrat<TGenoPheno>::_parameters._fplot.empty())
      {
	_fplotstream = new std::ofstream(eostrat<TGenoPheno>::_parameters._fplot);
	_fplotstream->precision(std::numeric_limits<double>::digits10);
      }
    auto mit=eostrat<TGenoPheno>::_parameters._stoppingcrit.begin();
    while(mit!=eostrat<TGenoPheno>::_parameters._stoppingcrit.end())
      {
	_stopcriteria.set_criteria_active((*mit).first,(*mit).second);
	++mit;
      }
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  CMAStrategy<TCovarianceUpdate,TGenoPheno>::CMAStrategy(FitFunc &func,
							 CMAParameters<TGenoPheno> &parameters,
							 const CMASolutions &cmasols)
    :ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >(func,parameters,cmasols)
  {
    eostrat<TGenoPheno>::_pfunc = _defaultPFunc;
    if (!parameters._full_fplot)
      eostrat<TGenoPheno>::_pffunc = _defaultFPFunc;
    else eostrat<TGenoPheno>::_pffunc = &fpfuncdef_full_impl<TCovarianceUpdate,TGenoPheno>;
    _esolver = Eigen::EigenMultivariateNormal<double>(false,eostrat<TGenoPheno>::_parameters._seed); // seeding the multivariate normal generator.
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
    if (!eostrat<TGenoPheno>::_parameters._sep && !eostrat<TGenoPheno>::_parameters._vd)
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
    else if (eostrat<TGenoPheno>::_parameters._sep)
      {
	_esolver.setMean(eostrat<TGenoPheno>::_solutions._xmean);
	_esolver.set_covar(eostrat<TGenoPheno>::_solutions._sepcov);
	_esolver.set_transform(eostrat<TGenoPheno>::_solutions._sepcov.cwiseSqrt());
      }
    else if (eostrat<TGenoPheno>::_parameters._vd)
      {
	_esolver.setMean(eostrat<TGenoPheno>::_solutions._xmean);
	_esolver.set_covar(eostrat<TGenoPheno>::_solutions._sepcov);
      }

    //debug
    //std::cout << "transform: " << _esolver._transform << std::endl;
    //debug
    
    // sample for multivariate normal distribution, produces one candidate per column.
    dMat pop;
    if (!eostrat<TGenoPheno>::_parameters._sep && !eostrat<TGenoPheno>::_parameters._vd)
      pop = _esolver.samples(eostrat<TGenoPheno>::_parameters._lambda,eostrat<TGenoPheno>::_solutions._sigma); // Eq (1).
    else if (eostrat<TGenoPheno>::_parameters._sep)
      pop = _esolver.samples_ind(eostrat<TGenoPheno>::_parameters._lambda,eostrat<TGenoPheno>::_solutions._sigma);
    else if (eostrat<TGenoPheno>::_parameters._vd)
      {
	pop = _esolver.samples_ind(eostrat<TGenoPheno>::_parameters._lambda);
	double normv = eostrat<TGenoPheno>::_solutions._v.squaredNorm();
	double fact = std::sqrt(1+normv)-1;
	dVec vbar = eostrat<TGenoPheno>::_solutions._v / std::sqrt(normv);

	pop += fact * vbar * (vbar.transpose() * pop);
	for (int i=0;i<pop.cols();i++)
	  {
	    pop.col(i) = eostrat<TGenoPheno>::_solutions._xmean + eostrat<TGenoPheno>::_solutions._sigma * eostrat<TGenoPheno>::_solutions._sepcov.cwiseProduct(pop.col(i));
	  }
      }
    
    // gradient if available.
    if (eostrat<TGenoPheno>::_parameters._with_gradient)
      {
	dVec grad_at_mean = eostrat<TGenoPheno>::gradf(eostrat<TGenoPheno>::_parameters._gp.pheno(eostrat<TGenoPheno>::_solutions._xmean));
	dVec gradgp_at_mean = eostrat<TGenoPheno>::gradgp(eostrat<TGenoPheno>::_solutions._xmean); // for geno / pheno transform.
	grad_at_mean = grad_at_mean.cwiseProduct(gradgp_at_mean);
	if (grad_at_mean != dVec::Zero(eostrat<TGenoPheno>::_parameters._dim))
	  {
	    dVec nx;
	    if (!eostrat<TGenoPheno>::_parameters._sep && !eostrat<TGenoPheno>::_parameters._vd)
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

    // tpa: fill up two first (or second in case of gradient) points with candidates usable for tpa computation
    if (eostrat<TGenoPheno>::_parameters._tpa == 2  && eostrat<TGenoPheno>::_niter > 0)
      {
	dVec mean_shift = eostrat<TGenoPheno>::_solutions._xmean - eostrat<TGenoPheno>::_solutions._xmean_prev;
	double mean_shift_norm = 1.0;
	if (!eostrat<TGenoPheno>::_parameters._sep && !eostrat<TGenoPheno>::_parameters._vd)
	  mean_shift_norm = (_esolver._eigenSolver.eigenvalues().cwiseSqrt().cwiseInverse().cwiseProduct(_esolver._eigenSolver.eigenvectors().transpose()*mean_shift)).norm() / eostrat<TGenoPheno>::_solutions._sigma;
	else mean_shift_norm = eostrat<TGenoPheno>::_solutions._sepcov.cwiseSqrt().cwiseInverse().cwiseProduct(mean_shift).norm() / eostrat<TGenoPheno>::_solutions._sigma;
	//std::cout << "mean_shift_norm=" << mean_shift_norm << " / sqrt(N)=" << std::sqrt(std::sqrt(eostrat<TGenoPheno>::_parameters._dim)) << std::endl;

	dMat rz = _esolver.samples_ind(1);
	double mfactor = rz.norm();
	dVec z = mfactor * (mean_shift / mean_shift_norm);
	eostrat<TGenoPheno>::_solutions._tpa_x1 = eostrat<TGenoPheno>::_solutions._xmean + z;
	eostrat<TGenoPheno>::_solutions._tpa_x2 = eostrat<TGenoPheno>::_solutions._xmean - z;
	
	// if gradient is in col 0, move tpa vectors to pos 1 & 2
	if (eostrat<TGenoPheno>::_parameters._with_gradient)
	  {
	    eostrat<TGenoPheno>::_solutions._tpa_p1 = 1;
	    eostrat<TGenoPheno>::_solutions._tpa_p2 = 2;
	  }
	pop.col(eostrat<TGenoPheno>::_solutions._tpa_p1) = eostrat<TGenoPheno>::_solutions._tpa_x1;
	pop.col(eostrat<TGenoPheno>::_solutions._tpa_p2) = eostrat<TGenoPheno>::_solutions._tpa_x2;
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

#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
#endif
    
    // sort candidates.
    if (!eostrat<TGenoPheno>::_parameters._uh)
      eostrat<TGenoPheno>::_solutions.sort_candidates();
    else eostrat<TGenoPheno>::uncertainty_handling();
    
    // call on tpa computation of s(t)
    if (eostrat<TGenoPheno>::_parameters._tpa == 2 && eostrat<TGenoPheno>::_niter > 0)
      eostrat<TGenoPheno>::tpa_update();

    // update function value history, as needed.
    eostrat<TGenoPheno>::_solutions.update_best_candidates();
    
    // CMA-ES update, depends on the selected 'flavor'.
    TCovarianceUpdate::update(eostrat<TGenoPheno>::_parameters,_esolver,eostrat<TGenoPheno>::_solutions);
    
    if (eostrat<TGenoPheno>::_parameters._uh)
      if (eostrat<TGenoPheno>::_solutions._suh > 0.0)
	eostrat<TGenoPheno>::_solutions._sigma *= eostrat<TGenoPheno>::_parameters._alphathuh;

    // other stuff.
    if (!eostrat<TGenoPheno>::_parameters._sep && !eostrat<TGenoPheno>::_parameters._vd)
      eostrat<TGenoPheno>::_solutions.update_eigenv(_esolver._eigenSolver.eigenvalues(),
						    _esolver._eigenSolver.eigenvectors());
    else eostrat<TGenoPheno>::_solutions.update_eigenv(eostrat<TGenoPheno>::_solutions._sepcov,
						       dMat::Constant(eostrat<TGenoPheno>::_parameters._dim,1,1.0));
#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
    eostrat<TGenoPheno>::_solutions._elapsed_tell = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  bool CMAStrategy<TCovarianceUpdate,TGenoPheno>::stop()
  {
    if (eostrat<TGenoPheno>::_solutions._run_status < 0) // an error occured, most likely out of memory at cov matrix creation.
      return true;
    
    if (eostrat<TGenoPheno>::_pfunc(eostrat<TGenoPheno>::_parameters,eostrat<TGenoPheno>::_solutions)) // progress function.
      return true; // end on progress function internal termination, possibly custom.
    
    if (!eostrat<TGenoPheno>::_parameters._fplot.empty())
      plot();
    
    if (eostrat<TGenoPheno>::_niter == 0)
      return false;

    if ((eostrat<TGenoPheno>::_solutions._run_status = _stopcriteria.stop(eostrat<TGenoPheno>::_parameters,eostrat<TGenoPheno>::_solutions)) != CONT)
      return true;
    else return false;
  }

  template <class TCovarianceUpdate, class TGenoPheno>
  int CMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize(const EvalFunc &evalf, const AskFunc &askf, const TellFunc &tellf)
  {
    //debug
    //DLOG(INFO) << "optimize()\n";
    //debug

    if (eostrat<TGenoPheno>::_initial_elitist 
	|| eostrat<TGenoPheno>::_parameters._initial_elitist
	|| eostrat<TGenoPheno>::_parameters._elitist
	|| eostrat<TGenoPheno>::_parameters._initial_fvalue)
      {
	eostrat<TGenoPheno>::_solutions._initial_candidate = Candidate(eostrat<TGenoPheno>::_func(eostrat<TGenoPheno>::_parameters._gp.pheno(eostrat<TGenoPheno>::_solutions._xmean).data(),eostrat<TGenoPheno>::_parameters._dim),
								       eostrat<TGenoPheno>::_solutions._xmean);
	eostrat<TGenoPheno>::_solutions._best_seen_candidate = eostrat<TGenoPheno>::_solutions._initial_candidate;
	this->update_fevals(1);
      }
    
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
    while(!stop())
      {
	dMat candidates = askf();
	evalf(candidates,eostrat<TGenoPheno>::_parameters._gp.pheno(candidates));
	tellf();
	eostrat<TGenoPheno>::inc_iter();
	std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	eostrat<TGenoPheno>::_solutions._elapsed_last_iter = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	tstart = std::chrono::system_clock::now();
      }
    if (eostrat<TGenoPheno>::_parameters._with_edm)
      eostrat<TGenoPheno>::edm();

    // test on final value wrt. to best candidate value and number of iterations in between.
    if (eostrat<TGenoPheno>::_parameters._initial_elitist_on_restart)
      {
	if (eostrat<TGenoPheno>::_parameters._initial_elitist_on_restart
	    && eostrat<TGenoPheno>::_solutions._best_seen_candidate.get_fvalue()
	    < eostrat<TGenoPheno>::_solutions.best_candidate().get_fvalue()
	    && eostrat<TGenoPheno>::_niter - eostrat<TGenoPheno>::_solutions._best_seen_iter >= 3) // elitist
	  {
	    LOG_IF(INFO,!eostrat<TGenoPheno>::_parameters._quiet) << "Starting elitist on restart: bfvalue=" << eostrat<TGenoPheno>::_solutions._best_seen_candidate.get_fvalue() << " / biter=" << eostrat<TGenoPheno>::_solutions._best_seen_iter << std::endl;
	    this->set_initial_elitist(true);

	    // reinit solution and re-optimize.
	    eostrat<TGenoPheno>::_parameters.set_x0(eostrat<TGenoPheno>::_solutions._best_seen_candidate.get_x_dvec_ref());
	    eostrat<TGenoPheno>::_solutions = CMASolutions(eostrat<TGenoPheno>::_parameters);
	    eostrat<TGenoPheno>::_solutions._nevals = eostrat<TGenoPheno>::_nevals;
	    eostrat<TGenoPheno>::_niter = 0;
	    optimize();
	  }
      }
    
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
  template class CMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class CMAStrategy<VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
