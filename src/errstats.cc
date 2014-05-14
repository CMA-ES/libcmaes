/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
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

#include "errstats.h"
#include <iostream>

namespace libcmaes
{
  
  template <class TGenoPheno>
  pli errstats<TGenoPheno>::profile_likelihood(FitFunc &func,
					       CMAParameters<TGenoPheno> &parameters,
					       CMASolutions &cmasol,
					       const int &k,
					       const bool &curve,
					       const int &samplesize,
					       const double &fup,
					       const double &delta)
  {
    dVec x = cmasol.best_candidate()._x;
    double minfvalue = cmasol.best_candidate()._fvalue;

    //debug
    //std::cout << "xk=" << x[k] << " / minfvalue=" << minfvalue << std::endl;
    //debug

    pli le(k,samplesize,parameters._dim,x,minfvalue);
    
    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,false,samplesize,fup,delta,curve); // positive direction
    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,true,samplesize,fup,delta,curve);  // negative direction

    le.setMinMax();
    cmasol._pls.insert(std::pair<int,pli>(k,le));
    return le;
  }

  template <class TGenoPheno>
  void errstats<TGenoPheno>::profile_likelihood_search(FitFunc &func,
						       CMAParameters<TGenoPheno> &parameters,
						       pli &le,
						       const CMASolutions &cmasol,
						       const int &k,
						       const bool &neg,
						       const int &samplesize,
						       const double &fup,
						       const double &delta,
						       const bool &curve)
  {
    int sign = neg ? -1 : 1;
    dVec x = cmasol.best_candidate()._x;
    double xk = x[k];
    double minfvalue = cmasol.best_candidate()._fvalue;
    CMASolutions citsol = cmasol;
    double dxk = sign * xk * 0.1;
    for (int i=0;i<samplesize;i++)
      {
	// get a new xk point.
	bool iterend = errstats<TGenoPheno>::take_linear_step(func,k,minfvalue,fup,curve,x,dxk);

	//debug
	//std::cout << "new xk point: " << x.transpose() << std::endl;
	//debug
	
	// minimize.
	citsol = errstats<TGenoPheno>::optimize_pk(func,parameters,citsol,k,x[k]);
	x = citsol.best_candidate()._x;
	minfvalue = citsol.best_candidate()._fvalue;
	
	// store points.
	le._fvaluem[samplesize+sign*(1+i)] = citsol.best_candidate()._fvalue;
	le._xm.row(samplesize+sign*(1+i)) = citsol.best_candidate()._x.transpose();

	if (!curve && iterend)
	  {
	    // pad and return.
	    for (int j=i+1;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = citsol.best_candidate()._fvalue;
		le._xm.row(samplesize+sign*(1+j)) = citsol.best_candidate()._x.transpose();
	      }
	    return;
	  }
      }
  }
						
  template <class TGenoPheno>
  bool errstats<TGenoPheno>::take_linear_step(FitFunc &func,
					      const int &k,
					      const double &minfvalue,
					      const double &fup,
					      const bool &curve,
					      dVec &x,
					      double &dxk)
  {
    static double fdiff_relative_increase = 0.1;
    double fdelta = 0.1 * fup;
    double threshold = minfvalue + fup;
    dVec xtmp = x;
    xtmp[k] += dxk;
    double fvalue = func(xtmp.data(),xtmp.size());
    double fdiff = fvalue - minfvalue;

    //debug
    //std::cout << "threshold=" << threshold << " / fdiff=" << fdiff << std::endl;
    //debug
    
    if (fdiff > threshold * fdiff_relative_increase) // decrease dxk
      {
	while((curve || fabs(fvalue-fup)>fdelta)
	      && fdiff > threshold * fdiff_relative_increase)
	  {
	    dxk /= 2.0;
	    xtmp[k] = x[k] + dxk;
	    fvalue = func(xtmp.data(),xtmp.size());
	    fdiff = fvalue - minfvalue;
	  }
      }
    else // increase dxk
      {
	while ((curve || fabs(fvalue-fup)>fdelta)
	       && fdiff < threshold * fdiff_relative_increase)
	  {
	    dxk *= 2.0;
	    xtmp[k] = x[k] + dxk;
	    fvalue = func(xtmp.data(),xtmp.size());
	    fdiff = fvalue - minfvalue;
	  }
	dxk /= 2.0;
      }
    x[k] += dxk; // set value.
    return (fabs(fvalue-fup) < fdelta);
  }

  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_pk(FitFunc &func,
						 CMAParameters<TGenoPheno> &parameters,
						 const CMASolutions &cmasol,
						 const int &k,
						 const double &vk)
  {
    CMASolutions ncmasol = cmasol;
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters._quiet = true; //TODO: option.
    nparameters.set_fixed_p(k,vk);
    return cmaes(func,nparameters);
  }
    
  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_reduced_pk(FitFunc &func,
							 CMAParameters<TGenoPheno> &parameters,
							 const CMASolutions &cmasol,
							 const int &k,
							 const double &vk)
  {
    // set re-arranged solution as starting point of the new optimization problem.
    CMASolutions ncmasol = cmasol;
    ncmasol.reset_as_fixed(k);

    // re-arrange function
    FitFunc nfunc = [func,k,vk](const double *x, const int N)
      {
	std::vector<double> nx;
	nx.reserve(N+1);
	for (int i=0;i<N+1;i++)
	  {
	    if (i < k)
	      nx.push_back(x[i]);
	    else if (i == k)
	      nx.push_back(vk);
	    else nx.push_back(x[i-1]);
	  }
	return func(&nx.front(),N+1);
      };
    
    // optimize and return result.
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters.reset_as_fixed(k);
    return cmaes(nfunc,nparameters);
  }

  template class errstats<GenoPheno<NoBoundStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy>>;
}
