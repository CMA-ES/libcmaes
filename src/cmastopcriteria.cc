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

#include "cmastopcriteria.h"
#include <math.h>
#include <iterator>
#include <glog/logging.h>
#include <limits>
#include <iostream>

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
  
  CMAStopCriteria::CMAStopCriteria()
    :_active(true)
  {
    StopCriteriaFunc autoMaxIter = [](const CMAParameters &cmap, const CMASolutions &cmas)
      {
	int thresh = static_cast<int>(100.0 + 50*pow(cmap._dim+3,2) / sqrt(cmap._lambda));
	if (cmas._niter >= thresh)
	  {
	    LOG_IF(INFO,!cmap._quiet) << "stopping criteria autoMaxIter => thresh=" << thresh << std::endl;
	    return AUTOMAXITER;
	  }
	return CONT;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(AUTOMAXITER,autoMaxIter));
    StopCriteriaFunc tolHistFun = [](const CMAParameters &cmap, const CMASolutions &cmas)
      {
	static double threshold = 1e-12;
	int histsize = static_cast<int>(cmas._best_candidates_hist.size());
	int histthresh = static_cast<int>(10+floor(30*cmap._dim/cmap._lambda));
	int histlength = std::min(histthresh,histsize);
	if (histlength < histthresh) // not enough data
	  return CONT;
	std::pair<double,double> frange(std::numeric_limits<double>::max(),std::numeric_limits<double>::min());
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(TOLHISTFUN,tolHistFun));
    StopCriteriaFunc equalFunVals = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(EQUALFUNVALS,equalFunVals));
    StopCriteriaFunc tolX = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
	  if (sqrt(cmas._cov(i,i))>tfactor)
	    return CONT;
	LOG_IF(INFO,!cmap._quiet) << "stopping criteria tolX\n";
	return TOLX;
      };
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(TOLX,tolX));
    StopCriteriaFunc tolUpSigma = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(TOLUPSIGMA,tolUpSigma));
    StopCriteriaFunc stagnation = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(STAGNATION,stagnation));

    StopCriteriaFunc conditionCov = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(CONDITIONCOV,conditionCov));
    StopCriteriaFunc noEffectAxis = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(NOEFFECTAXIS,noEffectAxis));
    StopCriteriaFunc noEffectCoor = [](const CMAParameters &cmap, const CMASolutions &cmas)
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
    _scriteria.insert(std::pair<int,StopCriteriaFunc>(NOEFFECTCOOR,noEffectCoor));
  }

  CMAStopCriteria::~CMAStopCriteria()
  {
  }

  int CMAStopCriteria::stop(const CMAParameters &cmap, const CMASolutions &cmas) const
  {
    if (!_active)
      return 0;
    int r = 0;
    for (auto imap : _scriteria)
      {
	if ((r=imap.second(cmap,cmas))!=0)
	  return r;
      }
    return 0;
  }
  
}
