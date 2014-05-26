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
#include "cmastopcriteria.h"
#include <cmath>
#include <iterator>
#include <glog/logging.h>
#include <limits>
#include <iostream>

#ifdef HAVE_DEBUG
#include <chrono>
#endif

namespace libcmaes
{

  // computes median of a vector.
  double median(std::vector<double> scores)
  {
    double median;
    size_t size = scores.size();
    std::sort(scores.begin(), scores.end());
    if (size  % 2 == 0)
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    else median = scores[size / 2];
    return median;
  }
  
  template <class TGenoPheno>
  CMAStopCriteria<TGenoPheno>::CMAStopCriteria()
    :_active(true)
  {
    StopCriteriaFunc<TGenoPheno> autoMaxIter = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	if (!cmap._has_max_iter) // this criteria is deactivated
	  return CONT;
	static int thresh = static_cast<int>(100.0 + 50*pow(cmap._dim+3,2) / sqrt(cmap._lambda));
	if (cmas._niter >= thresh)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria autoMaxIter => thresh=" << thresh << std::endl;
	    return AUTOMAXITER;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(AUTOMAXITER,autoMaxIter));
    StopCriteriaFunc<TGenoPheno> tolHistFun = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	static double threshold = std::max(cmap._ftolerance,1e-12); // set it once
	int histsize = static_cast<int>(cmas._best_candidates_hist.size());
	int histthresh = static_cast<int>(10+ceil(30*cmap._dim/cmap._lambda));
	int histlength = std::min(histthresh,histsize);
	if (histlength < histthresh) // not enough data
	  return CONT;
	std::pair<double,double> frange(std::numeric_limits<double>::max(),-std::numeric_limits<double>::max());
	for (int i=0;i<histlength;i++)
	  {
	    double val = cmas._best_candidates_hist.at(histsize-1-i)._fvalue;
	    frange.first = std::min(val,frange.first);
	    frange.second = std::max(val,frange.second);
	  }
	double rg = fabs(frange.second-frange.first);
	if (rg < threshold)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria tolHistFun => frange=" << rg << std::endl;
	    return TOLHISTFUN;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(TOLHISTFUN,tolHistFun));
    StopCriteriaFunc<TGenoPheno> equalFunVals = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	int histsize = static_cast<int>(cmas._best_candidates_hist.size());
	int histlength = std::min(cmap._dim,histsize);
	if (histlength < cmap._dim) // not enough data
	  return CONT;

	int c = 0;
	for (int i=0;i<histlength;i++)
	  {
	    if (cmas._best_candidates_hist.at(histsize-1-i)._fvalue
		== cmas._k_best_candidates_hist.at(histsize-1-i)._fvalue)
	      c++;
	  }
	if (c > histlength / 3.0)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria equalFunVals\n";
	    return EQUALFUNVALS;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(EQUALFUNVALS,equalFunVals));
    StopCriteriaFunc<TGenoPheno> tolX = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	static double tolx = 1e-12;
	double factor = cmas._sigma / cmap._sigma_init;
	double tfactor = tolx * factor;
	// test1: all components of pc . factor < tolx.
	for (int i=0;i<cmas._pc.rows();i++)
	  if (cmas._pc[i]>=tfactor)
	    return CONT;
	//test 2: all square root components of cov . factor < tolx.
	for (int i=0;i<cmas._cov.rows();i++)
	  if (sqrt(cmas._cov(i,i))>=tfactor)
	    return CONT;
	LOG_IF(INFO,!cmap._quiet) << "stopping criteria tolX\n";
	return TOLX;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(TOLX,tolX));
    StopCriteriaFunc<TGenoPheno> tolUpSigma = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	static double tolupsigma = 1e20;
	double factor = cmas._sigma / cmap._sigma_init;
	double rhs = tolupsigma * sqrt(cmas._max_eigenv);
	if (factor > rhs)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria tolUpSigma => factor=" << factor << " / max eigenv=" << cmas._max_eigenv << " / rhs=" << rhs << std::endl;
	    return TOLUPSIGMA;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(TOLUPSIGMA,tolUpSigma));
    StopCriteriaFunc<TGenoPheno> stagnation = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	if (cmas._bfvalues.size() < 20 || cmas._median_fvalues.size() < 20)
	  return CONT;
	double medianbv = median(cmas._bfvalues);
	std::vector<double> oldest_median_fvalues(20);
	std::copy_n(cmas._median_fvalues.begin(),20,oldest_median_fvalues.begin());
	double old_medianbv = median(oldest_median_fvalues);
	if (medianbv > old_medianbv)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria stagnation => oldmedianfvalue=" << old_medianbv << " / newmedianfvalue=" << medianbv << std::endl;
	    return STAGNATION;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(STAGNATION,stagnation));
    StopCriteriaFunc<TGenoPheno> conditionCov = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	static double bound = 1e14;
	double kappa = cmas._max_eigenv / cmas._min_eigenv;
	if (kappa > bound)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria conditionCov => min eigenv=" << cmas._min_eigenv << " / max eigenv=" << cmas._max_eigenv << " / kappa=" << kappa << std::endl;
	    return CONDITIONCOV;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(CONDITIONCOV,conditionCov));
    StopCriteriaFunc<TGenoPheno> noEffectAxis = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	double fact = 0.1*cmas._sigma;
	for (int i=0;i<cmap._dim;i++)
	  {
	    double ei = fact * sqrt(cmas._leigenvalues(i));
	    for (int j=0;j<cmap._dim;j++)
	      if (cmas._xmean[i] != cmas._xmean[i] + ei * cmas._leigenvectors(i,j))
		return CONT;
	  }
	LOG_IF(INFO,!cmap._quiet) << "stopping criteria NoEffectAxis\n";
	return NOEFFECTAXIS;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(NOEFFECTAXIS,noEffectAxis));
    StopCriteriaFunc<TGenoPheno> noEffectCoor = [](const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas)
      {
	double fact = 0.2*cmas._sigma;
	for (int i=0;i<cmap._dim;i++)
	  if (cmas._xmean[i] == fact * sqrt(cmas._cov(i,i)))
	    {
	      LOG_IF(INFO,!cmap._quiet) << "stopping criteria NoEffectCoor\n";
	      return NOEFFECTCOOR;
	    }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc<TGenoPheno>>(NOEFFECTCOOR,noEffectCoor));
  }

  template <class TGenoPheno>
  CMAStopCriteria<TGenoPheno>::~CMAStopCriteria()
  {
  }

  template <class TGenoPheno>
  int CMAStopCriteria<TGenoPheno>::stop(const CMAParameters<TGenoPheno> &cmap, const CMASolutions &cmas) const
  {
#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
#endif
    if (!_active)
      return 0;
    int r = 0;
    for (auto imap : _scriteria)
      {
	if ((r=imap.second(cmap,cmas))!=0)
	  {
#ifdef HAVE_DEBUG
	    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	    const_cast<CMASolutions&>(cmas)._elapsed_stop = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
	    return r;
	  }
      }
#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
    const_cast<CMASolutions&>(cmas)._elapsed_stop = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
    return 0;
  }

  template class CMAStopCriteria<GenoPheno<NoBoundStrategy>>;
  template class CMAStopCriteria<GenoPheno<pwqBoundStrategy>>;
  template class CMAStopCriteria<GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class CMAStopCriteria<GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
