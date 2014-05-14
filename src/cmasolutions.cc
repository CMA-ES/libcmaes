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
	_cov = dMat::Identity(p._dim,p._dim);
      }
    catch (std::bad_alloc &e)
      {
	_run_status = OPTI_ERR_OUTOFMEMORY;
	return;
      }
    if (p._x0min == p._x0max)
      {
	if (p._x0min == dVec::Constant(p._dim,std::numeric_limits<double>::min()))
	  _xmean = dVec::Random(p._dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
	else _xmean = p._x0min;
      }
    else
      {
	_xmean = 0.5*(dVec::Random(p._dim) + dVec::Constant(p._dim,1.0)); // scale to [0,1].
	_xmean = _xmean.cwiseProduct(p._x0max - p._x0min) + p._x0min; // scale to bounds.
      }
    if (static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init > 0.0)
      _sigma = static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init;
    else static_cast<CMAParameters<TGenoPheno>&>(p)._sigma_init = _sigma = 1.0/static_cast<double>(p._dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(p._dim);
    _pc = dVec::Zero(p._dim);
    _candidates.resize(p._lambda);
    _kcand = static_cast<int>(1.0+ceil(0.1+p._lambda/4.0));
  }

  CMASolutions::~CMASolutions()
  {
  }

  void CMASolutions::update_best_candidates()
  {
    _best_candidates_hist.push_back(_candidates.at(0)); // supposed candidates is sorted.
    _k_best_candidates_hist.push_back(_candidates.at(_kcand));

    _bfvalues.push_back(_candidates.at(0)._fvalue);
    if (_bfvalues.size() > 20)
      _bfvalues.erase(_bfvalues.begin());

    // get median of candidate's scores, used in termination criteria (stagnation).
    double median = 0.0;
    size_t csize = _candidates.size();
    if (csize % 2 == 0)
      median = (_candidates[csize/2-1]._fvalue + _candidates[csize/2]._fvalue)/2.0;
    else median = _candidates[csize/2]._fvalue;
    _median_fvalues.push_back(median);
    if (_median_fvalues.size() > static_cast<size_t>(ceil(0.2*_niter+120+30*_xmean.size()/static_cast<double>(_candidates.size()))))
      _median_fvalues.erase(_median_fvalues.begin());
    
    //debug
    /*std::cerr << "ordered candidates:\n";
    for (size_t i=0;i<_candidates.size();i++)
      {
	std::cerr << _candidates.at(i)._fvalue << " / " << _candidates.at(i)._x.transpose() << std::endl;
	}*/
    //debug
  }

  void CMASolutions::update_eigenv(const dVec &eigenvalues,
				   const dMat &eigenvectors)
  {
    _max_eigenv = eigenvalues.maxCoeff();
    _min_eigenv = eigenvalues.minCoeff();
    _leigenvalues = eigenvalues;
    _leigenvectors = eigenvectors;
  }

  void CMASolutions::reset_as_fixed(const int &k)
  {
    removeRow(_cov,k);
    removeColumn(_cov,k);
    removeRow(_csqinv,k);
    removeColumn(_csqinv,k);
    removeElement(_xmean,k);
    removeElement(_psigma,k);
    removeElement(_pc,k);
    for (size_t i=0;i<_candidates.size();i++)
      removeElement(_candidates.at(i)._x,k);
    _best_candidates_hist.clear();
    removeElement(_leigenvalues,k);
    removeRow(_leigenvectors,k);
    removeColumn(_leigenvectors,k);
    _niter = 0;
    _nevals = 0;
    _k_best_candidates_hist.clear();
    _bfvalues.clear();
    _median_fvalues.clear();
    _run_status = 0;
    _elapsed_time = _elapsed_last_iter = 0;
#ifdef HAVE_DEBUG
    _elapsed_eval = _elapsed_ask = _elapsed_tell = _elapsed_stop = 0;
#endif
  }
  
  std::ostream& CMASolutions::print(std::ostream &out,
				    const int &verb_level) const
  {
    out << "best solution => f-value=" << best_candidate()._fvalue << " / sigma=" << _sigma << " / iter=" << _niter << " / elaps=" << _elapsed_time << "ms" << " / x=" << best_candidate()._x.transpose(); //TODO: print pheno(x).
    if (verb_level)
      {
	out << "\ncovdiag=" << _cov.diagonal().transpose() << std::endl;
	out << "psigma=" << _psigma.transpose() << std::endl;
	out << "pc=" << _pc.transpose() << std::endl;
      }
    if (!_pls.empty())
      {
	out << "\nconfidence intervals:\n";
	for (auto it=_pls.begin();it!=_pls.end();++it)
	  {
	    out << "dim " << (*it).first << " in [" << (*it).second._min << "," << (*it).second._max << "]\n";
	  }
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
}
