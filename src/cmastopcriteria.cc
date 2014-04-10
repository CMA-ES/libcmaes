
#include "cmastopcriteria.h"
#include <math.h>
#include <iterator>
#include <glog/logging.h>
#include <iostream>

namespace libcmaes
{

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
	std::pair<double,double> frange(DBL_MAX,DBL_MIN);
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
	//TODO: requires median func values to be stored.
	
	return CONT;
      };

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
	
	return CONT;
      };
    StopCriteriaFunc noEffectCoor = [](const CMAParameters &cmap, const CMASolutions &cmas)
      {
	return CONT;
      };
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
	if ((r=imap.second(cmap,cmas))>0)
	  return r;
      }
    return 0;
  }
  
}
