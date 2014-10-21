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

#include "cmasolutions.h"
#include "opti_err.h"
#include "eigenmvn.h"
#include <limits>
#include <iostream>

namespace libcmaes
{
  template <class TGenoPheno>
  CMASolutions::CMASolutions(Parameters<TGenoPheno> &p)
    :_hsig(1),_max_eigenv(0.0),_min_eigenv(0.0),_niter(0),_nevals(0),_kcand(1),_eigeniter(0),_updated_eigen(true),_run_status(0),_elapsed_time(0)
  {
    try
      {
	if (!static_cast<CMAParameters<TGenoPheno>&>(p)._sep && !static_cast<CMAParameters<TGenoPheno>&>(p)._vd)
	  _cov = dMat::Identity(p._dim,p._dim);
	else _sepcov = dMat::Constant(p._dim,1,1.0);
      }
    catch (std::bad_alloc &e)
      {
	_run_status = OPTI_ERR_OUTOFMEMORY;
	return;
      }
    if (p._x0min == p._x0max)
      {
	if (p._x0min == dVec::Constant(p._dim,-std::numeric_limits<double>::max()))
	  _xmean = dVec::Random(p._dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
	else _xmean = p._x0min;
      }
    else
      {
	_xmean = 0.5*(dVec::Random(p._dim) + dVec::Constant(p._dim,1.0)); // scale to [0,1].
	_xmean = _xmean.cwiseProduct(p._x0max - p._x0min) + p._x0min; // scale to bounds.
      }
    if (!p._fixed_p.empty())
      {
	auto fpmit = p._fixed_p.begin();
	while (fpmit!=p._fixed_p.end())
	  {
	    _xmean((*fpmit).first) = (*fpmit).second;
	    ++fpmit;
	  }
      }
    // if scaling, need to apply to xmean.
    if (!p._gp._scalingstrategy._id)
      p._gp._scalingstrategy.scale_to_internal(_xmean,_xmean);
    if (static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init > 0.0)
      _sigma = static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init;
    else static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init = _sigma = 1.0/static_cast<double>(p._dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(p._dim);
    _pc = dVec::Zero(p._dim);
    _candidates.resize(p._lambda);
    _kcand = std::min(p._lambda-1,static_cast<int>(1.0+ceil(0.1+p._lambda/4.0)));
    _max_hist = (p._max_hist > 0) ? p._max_hist : static_cast<int>(10+ceil(30*p._dim/p._lambda));
    
    if (static_cast<CMAParameters<TGenoPheno>&>(p)._vd)
      {
	EigenMultivariateNormal<double> esolver(false,static_cast<uint64_t>(time(nullptr)));
	esolver.set_covar(_sepcov);
	_v = esolver.samples_ind(1) / std::sqrt(p._dim);
      }
  }

  CMASolutions::~CMASolutions()
  {
  }

  void CMASolutions::update_best_candidates()
  {
    _best_candidates_hist.push_back(_candidates.at(0)); // supposed candidates is sorted.
    _k_best_candidates_hist.push_back(_candidates.at(_kcand));
    if ((int)_best_candidates_hist.size() > _max_hist)
      {
	_best_candidates_hist.erase(_best_candidates_hist.begin());
	_k_best_candidates_hist.erase(_k_best_candidates_hist.begin());
      }
    
    _bfvalues.push_back(_candidates.at(0).get_fvalue());
    if (_bfvalues.size() > 20)
      _bfvalues.erase(_bfvalues.begin());

    // get median of candidate's scores, used in termination criteria (stagnation).
    double median = 0.0;
    size_t csize = _candidates.size();
    if (csize % 2 == 0)
      median = (_candidates[csize/2-1].get_fvalue() + _candidates[csize/2].get_fvalue())/2.0;
    else median = _candidates[csize/2].get_fvalue();
    _median_fvalues.push_back(median);
    if (_median_fvalues.size() > static_cast<size_t>(ceil(0.2*_niter+120+30*_xmean.size()/static_cast<double>(_candidates.size()))))
      _median_fvalues.erase(_median_fvalues.begin());
  }

  void CMASolutions::update_eigenv(const dVec &eigenvalues,
				   const dMat &eigenvectors)
  {
    _max_eigenv = eigenvalues.maxCoeff();
    _min_eigenv = eigenvalues.minCoeff();
    _leigenvalues = eigenvalues;
    _leigenvectors = eigenvectors;
  }

  template <class TGenoPheno>
  std::ostream& CMASolutions::print(std::ostream &out,
				    const int &verb_level,
				    const TGenoPheno &gp) const
  {
    if (_candidates.empty())
      {
	out << "empth solution set\n";
	return out;
      }
    out << "best solution => f-value=" << best_candidate().get_fvalue() << " / fevals=" << _nevals << " / sigma=" << _sigma << " / iter=" << _niter << " / elaps=" << _elapsed_time << "ms" << " / x=" << gp.pheno(best_candidate().get_x_dvec()).transpose();
    if (verb_level)
      {
	out << "\ncovdiag=" << _cov.diagonal().transpose() << std::endl;
	out << "psigma=" << _psigma.transpose() << std::endl;
	out << "pc=" << _pc.transpose() << std::endl;
      }
    return out;
  }

  std::ostream& operator<<(std::ostream &out, const CMASolutions &cmas)
  {
    cmas.print(out,0);
        return out;
  }
  
  template CMASolutions::CMASolutions(Parameters<GenoPheno<NoBoundStrategy>>&);
  template CMASolutions::CMASolutions(Parameters<GenoPheno<pwqBoundStrategy>>&);
  template CMASolutions::CMASolutions(Parameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>&);
  template CMASolutions::CMASolutions(Parameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>&);

  template std::ostream& CMASolutions::print(std::ostream&,const int&,const GenoPheno<NoBoundStrategy>&) const;
  template std::ostream& CMASolutions::print(std::ostream&,const int&,const GenoPheno<pwqBoundStrategy>&) const;
  template std::ostream& CMASolutions::print(std::ostream&,const int&,const GenoPheno<NoBoundStrategy,linScalingStrategy>&) const;
  template std::ostream& CMASolutions::print(std::ostream&,const int&,const GenoPheno<pwqBoundStrategy,linScalingStrategy>&) const;
}
