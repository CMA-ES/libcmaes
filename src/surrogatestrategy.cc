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

#include "surrogatestrategy.h"
#include "ipopcmastrategy.h"
#include "bipopcmastrategy.h"
#include <unordered_set>

namespace libcmaes
{

  template <class TGenoPheno> using eostrat = ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >;
  
  /*- SurrogateStrategy -*/
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::SurrogateStrategy(FitFunc &func, CMAParameters<TGenoPheno> &parameters)
    :TStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters)
  {
    _l = std::floor(30*std::sqrt(eostrat<TGenoPheno>::_parameters.dim()));
    eostrat<TGenoPheno>::_pfunc = [this](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols)
      {
	LOG_IF(INFO,!cmaparams.quiet()) << "iter=" << cmasols.niter() << " / evals=" << cmasols.fevals() << " / f-value=" << cmasols.best_candidate().get_fvalue() <<  " / sigma=" << cmasols.sigma() << " / trainerr=" << _train_err << " / testerr=" << _test_err << " / smtesterr=" << _smooth_test_err << std::endl;
	return 0;
      };
    eostrat<TGenoPheno>::_pffunc = [this](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
      {
	std::string sep = " ";
	fplotstream << fabs(cmasols.best_candidate().get_fvalue()) << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << sqrt(cmasols.max_eigenv()/cmasols.min_eigenv()) << sep;
	if (!cmasols.eigenvalues().size())
	  fplotstream << dVec::Zero(cmaparams.dim()).transpose() << sep;
	else fplotstream << cmasols.eigenvalues().transpose() << sep;
	if (!cmaparams.is_sep() && !cmaparams.is_vd())
	  fplotstream << cmasols.cov().sqrt().diagonal().transpose() << sep; // max deviation in all main axes
	else if (cmaparams.is_sep())
	  fplotstream << cmasols.sepcov().cwiseSqrt().transpose() << sep;
	else if (cmaparams.is_vd())
	  fplotstream << cmasols.sepcov().transpose() << sep;
	fplotstream << cmaparams.get_gp().pheno(cmasols.xmean()).transpose();
	fplotstream << sep << cmasols.elapsed_last_iter();
	fplotstream << sep << _train_err << sep << _test_err << sep << _smooth_test_err;
	fplotstream << std::endl;
	return 0;
      };
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::~SurrogateStrategy()
  {
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::add_to_training_set(const Candidate &c)
  {
    _tset.push_back(c);
    if ((int)_tset.size() > _l)
      _tset.erase(_tset.begin());
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  double SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::compute_error(const std::vector<Candidate> &test_set,
									const dMat &cov)
  {
    double test_set_size = test_set.size();
    double factor = 2.0 / static_cast<double>((test_set_size*(test_set_size-1)));
    std::vector<Candidate> ctest_set = test_set;

    std::sort(ctest_set.begin(),ctest_set.end(),
	      [](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() > c2.get_fvalue();});

    this->predict(ctest_set,cov);
    
    double err = 0.0;
    for (size_t i=0;i<ctest_set.size();i++)
      for (size_t j=i+1;j<ctest_set.size();j++)
	{
	  err += (ctest_set.at(i).get_fvalue() >= ctest_set.at(j).get_fvalue()) ? 0.0 : 1.0;
	}
    err *= factor;
    return err;
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::set_test_error(const double &err)
  {
    _test_err = err;
    _smooth_test_err = (1.0-_beta_err)*_smooth_test_err + _beta_err*err;
  }
  
  /*- SimpleSurrogateStrategy -*/
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::SimpleSurrogateStrategy(FitFunc &func,
										 CMAParameters<TGenoPheno> &parameters)
    :SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>(func,parameters)
  {
    eostrat<TGenoPheno>::_pfunc = [this](const CMAParameters<TGenoPheno> &cmaparams, const CMASolutions &cmasols)
      {
	LOG_IF(INFO,!cmaparams.quiet()) << "iter=" << cmasols.niter() << " / evals=" << cmasols.fevals() << " / f-value=" << cmasols.best_candidate().get_fvalue() <<  " / sigma=" << cmasols.sigma() << " / trainerr=" << this->_train_err << " / testerr=" << this->_test_err << " / smtesterr=" << this->_smooth_test_err << " / slifel=" << this->_nsteps << std::endl;
	return 0;
      };
    this->_stopcriteria.set_criteria_active(STAGNATION,false); // deactivate stagnation check due to the presence of ranks as median objective function values
    this->_stopcriteria.set_criteria_active(AUTOMAXITER,false);
  }
  
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::~SimpleSurrogateStrategy()
  {
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::eval(const dMat &candidates,
								   const dMat &phenocandidates)
  {
    if (!this->_exploit || this->_niter % this->_nsteps == 0 || (int)this->_tset.size() < this->_l)
      {
	// reactivate fvalue-based stopping criteria
	this->_stopcriteria.set_criteria_active(FTARGET,true);
	
	// compute test error if needed.
	if (this->_niter != 0 && (int)this->_tset.size() >= this->_l)
	  {
	    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
	      this->set_test_error(this->compute_error(eostrat<TGenoPheno>::_solutions._candidates,
						       eostrat<TGenoPheno>::_solutions._csqinv));
	    else this->set_test_error(this->compute_error(eostrat<TGenoPheno>::_solutions._candidates,
							  eostrat<TGenoPheno>::_solutions._sepcsqinv));
	    if (this->_auto_nsteps)
	      this->_nsteps = compute_lifelength();
	  }
	
	// use original objective function and collect points.
	eostrat<TGenoPheno>::eval(candidates,phenocandidates);
	for (int i=0;i<eostrat<TGenoPheno>::_solutions.size();i++)
	  this->add_to_training_set(eostrat<TGenoPheno>::_solutions.candidates().at(i)); // XXX: not very efficient update when test size already filled up.
      }
    else
      {
	// deactivate value-based stopping criteria
	this->_stopcriteria.set_criteria_active(FTARGET,false);
	
	// exploit surrogate
	for (int r=0;r<candidates.cols();r++)
	  {
	    eostrat<TGenoPheno>::_solutions.get_candidate(r).set_x(candidates.col(r));
	    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
	      this->predict(eostrat<TGenoPheno>::_solutions.candidates(),
			    eostrat<TGenoPheno>::_solutions.csqinv());
	    else this->predict(eostrat<TGenoPheno>::_solutions.candidates(),
			       eostrat<TGenoPheno>::_solutions.sepcsqinv());
	  }
      }
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::tell()
  {
    TStrategy<TCovarianceUpdate,TGenoPheno>::tell();

    // train surrogate as required.
    if (do_train())
      {
	this->train(this->_tset,
		    eostrat<TGenoPheno>::_solutions._csqinv);
      }
  }
  
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  int SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::optimize()
  {
    return TStrategy<TCovarianceUpdate,TGenoPheno>::optimize(std::bind(&SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::eval,this,std::placeholders::_1,std::placeholders::_2),
							     std::bind(&CMAStrategy<TCovarianceUpdate,TGenoPheno>::ask,this),
							     std::bind(&SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::tell,this));
  }
  
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  int SimpleSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::compute_lifelength()
  {
    int nsteps = std::max(1,static_cast<int>(std::floor((_terr - this->_smooth_test_err)/_terr * _nmax)));
    return nsteps+1;
  }

  /*- ACMSurrogateStrategy -*/
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::ACMSurrogateStrategy(FitFunc &func, CMAParameters<TGenoPheno> &parameters)
    :SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>(func,parameters)
  {
    _lambdaprime = std::floor(eostrat<TGenoPheno>::_parameters.lambda()/3.0);
    _prelambda = 50 * eostrat<TGenoPheno>::_parameters.lambda();
    init_rd();
    this->_stopcriteria.set_criteria_active(STAGNATION,false); // deactivate stagnation check due to the presence of ranks as median objective function values
  }
  
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::~ACMSurrogateStrategy()
  {
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::init_rd()
  {
    _gen0 = std::mt19937(_rd());
    _gen1 = std::mt19937(_rd());
    _norm_sel0 = std::normal_distribution<double>(0.0,_theta_sel0*_theta_sel0);
    _norm_sel1 = std::normal_distribution<double>(0.0,_theta_sel1*_theta_sel1);
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  dMat ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::ask()
  {
    // when starting or restarting, make sure the training set is reset.
    if (this->_niter == 0)
      {
	if (!this->_tset.empty()) // this is a restart
	  _prelambda = 50 * eostrat<TGenoPheno>::_parameters._lambda; // rescale pre-lambda, in case this is the first step of a restart with different population, also this is overriding the default upon restart
	this->reset_training_set();
      }
	
    if (this->_exploit && (int)this->_tset.size() >= this->_l)
      {
	double lambda = eostrat<TGenoPheno>::_parameters._lambda; // XXX: hacky.
	eostrat<TGenoPheno>::_parameters._lambda = _prelambda;
	dMat pop = TStrategy<TCovarianceUpdate,TGenoPheno>::ask();
	eostrat<TGenoPheno>::_parameters._lambda = lambda;
	return pop;
      }
    else return TStrategy<TCovarianceUpdate,TGenoPheno>::ask();
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::eval(const dMat &candidates,
									  const dMat &phenocandidates)
  {
    // use pre selection eval only if surrogate is ready.
    if (this->_exploit && (int)this->_tset.size() >= this->_l)
      {
	pre_selection_eval(candidates);
      }
    else
      {
	// compute test error if needed.
	if (this->_niter != 0 && (int)this->_tset.size() >= this->_l)
	  {
	    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
	      this->set_test_error(this->compute_error(eostrat<TGenoPheno>::_solutions.candidates(),
						       eostrat<TGenoPheno>::_solutions._csqinv));
	    else this->set_test_error(this->compute_error(eostrat<TGenoPheno>::_solutions.candidates(),
							  eostrat<TGenoPheno>::_solutions._sepcsqinv));
	  }
	
	// use original objective function and collect points.
	eostrat<TGenoPheno>::eval(candidates,phenocandidates);
	for (int i=0;i<eostrat<TGenoPheno>::_solutions.size();i++)
	  this->add_to_training_set(eostrat<TGenoPheno>::_solutions.candidates().at(i)); // XXX: not very efficient update when test size already filled up.
      }
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::tell()
  {
    if (!this->_exploit || (int)this->_tset.size() < this->_l)
      eostrat<TGenoPheno>::_solutions.sort_candidates();
    
    // update function value history, as needed.
    eostrat<TGenoPheno>::_solutions.update_best_candidates();
        
    // CMA-ES update, depends on the selected 'flavor'.
    TCovarianceUpdate::update(eostrat<TGenoPheno>::_parameters,this->_esolver,eostrat<TGenoPheno>::_solutions);
    
    // other stuff.
    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
      eostrat<TGenoPheno>::_solutions.update_eigenv(this->_esolver._eigenSolver.eigenvalues(),
						    this->_esolver._eigenSolver.eigenvectors());
    else eostrat<TGenoPheno>::_solutions.update_eigenv(eostrat<TGenoPheno>::_solutions._sepcov,
						       dMat::Constant(eostrat<TGenoPheno>::_parameters._dim,1,1.0));
    
    // train surrogate as required.
    if (do_train())
      {
	if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
	  {
	    this->train(this->_tset,eostrat<TGenoPheno>::_solutions._csqinv);
	    this->set_train_error(this->compute_error(this->_tset,eostrat<TGenoPheno>::_solutions._csqinv));
	  }
	else 
	  {
	    this->train(this->_tset,eostrat<TGenoPheno>::_solutions._sepcsqinv);
	    this->set_train_error(this->compute_error(this->_tset,eostrat<TGenoPheno>::_solutions._sepcsqinv));
	  }
      }
  }

  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  int ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::optimize()
  {
    return TStrategy<TCovarianceUpdate,TGenoPheno>::optimize(std::bind(&ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::eval,this,std::placeholders::_1,std::placeholders::_2),
							     std::bind(&ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::ask,this),
							     std::bind(&ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::tell,this));
  }
  
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate, class TGenoPheno>
  void ACMSurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::pre_selection_eval(const dMat &candidates)
  {
    // - rank all candidates according to surrogate.
    eostrat<TGenoPheno>::_solutions._candidates.clear();
    for (int r=0;r<candidates.cols();r++)
      {
	eostrat<TGenoPheno>::_solutions._candidates.push_back(Candidate(0.0,candidates.col(r)));
      }
    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
      this->predict(eostrat<TGenoPheno>::_solutions._candidates,eostrat<TGenoPheno>::_solutions._csqinv);
    else this->predict(eostrat<TGenoPheno>::_solutions._candidates,eostrat<TGenoPheno>::_solutions._sepcsqinv);
    eostrat<TGenoPheno>::_solutions.sort_candidates();

    /*std::vector<Candidate> &vc = eostrat<TGenoPheno>::_solutions.candidates();
    for (size_t i=0;i<vc.size();i++)
      {
	std::cout << vc.at(i).get_fvalue() << " / x=" << vc.at(i).get_x_dvec().transpose() << std::endl;
	}*/

    // - draw 'a'<lambda_pre samples according to lambda_pre*N(0,theta_sel0^2) and retain each sample from initial population, with rank r < floor(a)
    std::vector<Candidate> ncandidates;
    std::unordered_set<int> uh;
    
    // keep the estimated best candidate.
    ncandidates.push_back(eostrat<TGenoPheno>::_solutions._candidates.at(0));
    ncandidates.at(0).set_rank(0);
    uh.insert(0);
    
    std::unordered_set<int>::const_iterator uhit;
    while((int)ncandidates.size() < eostrat<TGenoPheno>::_parameters._lambda)
      {
	double da = _prelambda*std::fabs(_norm_sel0(_gen0));
	int a = std::floor(da);
	if (a < (int)eostrat<TGenoPheno>::_solutions._candidates.size() && (uhit=uh.find(a))==uh.end())
	  {
	    uh.insert(a);
	    eostrat<TGenoPheno>::_solutions._candidates.at(a).set_rank(a);
	    ncandidates.push_back(eostrat<TGenoPheno>::_solutions._candidates.at(a));
	  }
      }
    
    // - draw 'a'<lambda samples according to lambda*N(0,theta_sel1^2) and retain each sample from the population from previous step, with rank r < floor(a)
    std::vector<Candidate> test_set;
    std::sort(ncandidates.begin(),ncandidates.end(),
	      [](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() < c2.get_fvalue();});
    ncandidates.at(0).set_fvalue(this->_func(eostrat<TGenoPheno>::_parameters._gp.pheno(ncandidates.at(0).get_x_dvec()).data(),ncandidates.at(0).get_x_size()));
    this->update_fevals(1);
    test_set.push_back(ncandidates.at(0));
    this->add_to_training_set(ncandidates.at(0));
    int count = 1;
    uh.clear();
    uh.insert(0);
    while(count < _lambdaprime)
      {
	// XXX: do we need to drop the samples with a > lambda and if a == 0 ? weird bias on sampling...
	double da = std::fabs(_norm_sel1(_gen1));
	int a = std::floor(eostrat<TGenoPheno>::_parameters._lambda*da);
	if (a < (int)ncandidates.size() && (uhit=uh.find(a))==uh.end())
	  {
	    uh.insert(a);
	    double fvalue = this->_func(eostrat<TGenoPheno>::_parameters._gp.pheno(ncandidates.at(a).get_x_dvec()).data(),ncandidates.at(a).get_x_size());
	    ncandidates.at(a).set_fvalue(fvalue);
	    test_set.push_back(ncandidates.at(a));
	    this->add_to_training_set(ncandidates.at(a));
	    this->update_fevals(1);
	    ++count;
	  }
      }
    for (size_t i=1;i<ncandidates.size();i++)
      if ((uhit=uh.find(i))==uh.end())
	ncandidates.at(i).set_fvalue(std::numeric_limits<double>::max());
    std::sort(ncandidates.begin(),ncandidates.end(),
	      [](Candidate const &c1, Candidate const &c2)
	      {
		if (c1.get_fvalue() < c2.get_fvalue())
		  return true;
		if (c1.get_fvalue() == c2.get_fvalue())
		  return c1.get_rank() < c2.get_rank();
		return false;
	      });
    
    // test error.
    if (!eostrat<TGenoPheno>::_parameters.is_sep() && !eostrat<TGenoPheno>::_parameters.is_vd())
      this->set_test_error(this->compute_error(test_set,eostrat<TGenoPheno>::_solutions._csqinv));
    else this->set_test_error(this->compute_error(test_set,eostrat<TGenoPheno>::_solutions._sepcsqinv));
    
    // set candidate set.
    eostrat<TGenoPheno>::_solutions._candidates = ncandidates;
  }
  
  template class SimpleSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;

  template class SimpleSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;

  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class SimpleSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  
  template class ACMSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<CMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;

  template class ACMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<IPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;

  template class ACMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,CovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,ACovarianceUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
  template class ACMSurrogateStrategy<BIPOPCMAStrategy,VDCMAUpdate,GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
