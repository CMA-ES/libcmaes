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

#ifndef GENOPHENO_H
#define GENOPHENO_H

#include "noboundstrategy.h"
#include "pwq_bound_strategy.h"

namespace libcmaes
{
  typedef std::function<double (const double*, double*)> TransFunc;
  
  template <class TBoundStrategy=NoBoundStrategy>
    class GenoPheno
    {
    public:
    GenoPheno() {};
    GenoPheno(const double *lbounds, const double *ubounds, const int &dim)
    :_boundstrategy(lbounds,ubounds,dim)
    {};
    ~GenoPheno() {};

    dMat pheno(const dMat &candidates)
    {
      //TODO: apply custom pheno function.
      
      dMat ycandidates = dMat(candidates.rows(),candidates.cols());
#pragma omp parallel for if (candidates.cols() >= 100)
      for (int i=0;i<candidates.cols();i++)
	{
	  // apply bounds.
	  dVec xcoli = candidates.col(i);
	  dVec ycoli = ycandidates.col(i); // copy required as Eigen complains on function signature otherwise.
	  _boundstrategy.to_f_representation(xcoli,ycoli);
	  ycandidates.col(i) = ycoli;
	}
      return ycandidates;
    }

    dVec pheno(const dVec &candidate)
    {
      //TODO: apply custom pheno function.
      
      dVec phen = dVec::Zero(candidate.rows());
      _boundstrategy.to_f_representation(candidate,phen);
      return phen;
    }
    
    dVec geno(const dVec &candidate)
    {
      // reverse bounds.
      dVec gen = dVec::Zero(candidate.rows());
      _boundstrategy.to_internal_representation(gen,candidate);

      //TODO: apply custom geno function.

      return gen;
    }
    
    TBoundStrategy _boundstrategy;
  };

  // specialization when no bound strategy applies.
  template<> inline dMat GenoPheno<NoBoundStrategy>::pheno(const dMat &candidates)
    {
      return candidates;
    }
  template<> inline dVec GenoPheno<NoBoundStrategy>::pheno(const dVec &candidate)
    {
      return candidate;
    }
  template<> inline dVec GenoPheno<NoBoundStrategy>::geno(const dVec &candidate)
    {
      return candidate;
    }
}

#endif
