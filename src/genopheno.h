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

#ifndef GENOPHENO_H
#define GENOPHENO_H

#include "noboundstrategy.h"
#include "pwq_bound_strategy.h"
#include "scaling.h"

namespace libcmaes
{
  typedef std::function<void (const double*, double*, const int&)> TransFunc;

  /*TransFunc id_genof = [](const double *ext, const double *in)
    {
      in = ext;
    };

  TransFunc id_phenof = [](const double *in, const double *ext)
    {
      ext = in;
      };*/
  
  template <class TBoundStrategy=NoBoundStrategy,class TScalingStrategy=NoScalingStrategy>
    class GenoPheno
    {
    public:
    GenoPheno()
    :_id(true)
    {};

    GenoPheno(TransFunc &genof, TransFunc &phenof)
    :_genof(genof),_phenof(phenof),_id(false)
    {};
    
    GenoPheno(const double *lbounds, const double *ubounds, const int &dim)
    :_boundstrategy(lbounds,ubounds,dim),_id(true),_scalingstrategy(lbounds,ubounds,dim)
    {
      if (!_scalingstrategy._id)
	_boundstrategy = TBoundStrategy(lbounds,ubounds,dim);
      else
	{
	  std::vector<double> lb(dim,_scalingstrategy._intmin);
	  std::vector<double> ub(dim,_scalingstrategy._intmax);
	  _boundstrategy = TBoundStrategy(&lb.front(),&ub.front(),dim);
	}
    };

    GenoPheno(TransFunc &genof, TransFunc &phenof,
	      const double *lbounds, const double *ubounds, const int &dim)
    :_boundstrategy(lbounds,ubounds),_genof(genof),_phenof(phenof),_id(false),_scalingstrategy(lbounds,ubounds,dim)
    {
      if (!_scalingstrategy._id)
	_boundstrategy = TBoundStrategy(lbounds,ubounds,dim);
      else
	{
	  std::vector<double> lb(dim,_scalingstrategy._intmin);
	  std::vector<double> ub(dim,_scalingstrategy._intmax);
	  _boundstrategy = TBoundStrategy(&lb.front(),&ub.front(),dim);
	}
    };
    
    ~GenoPheno() {};

    private:
    dMat pheno_candidates(const dMat &candidates)
    {
      dMat ncandidates;
      if (!_id)
	{
	  ncandidates = dMat(candidates.rows(),candidates.cols());
#pragma omp parallel for if (candidates.cols() >= 100)
	  for (int i=0;i<candidates.cols();i++)
	    {
	      dVec ext = dVec(candidates.rows());
	      _phenof(candidates.col(i).data(),ext.data(),candidates.rows());
	      ncandidates.col(i) = ext;
	    }
	  return ncandidates;
	}
      return candidates;
    }
    
    public:
    dMat pheno(const dMat &candidates)
    {
      // apply custom pheno function.
      dMat ncandidates = pheno_candidates(candidates);
      
      // apply scaling.
      dMat ycandidates = dMat(ncandidates.rows(),ncandidates.cols());
      if (!_scalingstrategy._id)
	{
#pragma omp parallel for if (candidates.cols() >= 100)
	  for (int i=0;i<ncandidates.cols();i++)
	    {
	      dVec xcoli = ncandidates.col(i);
	      dVec ycoli;
	      _scalingstrategy.scale_to_f(xcoli,ycoli);
	      ycandidates.col(i) = ycoli;
	    }
	}

      // apply bounds.
#pragma omp parallel for if (candidates.cols() >= 100)
      for (int i=0;i<candidates.cols();i++)
	{
	  dVec xcoli;
	  if (_scalingstrategy._id)
	    xcoli = ncandidates.col(i);
	  else xcoli = ycandidates.col(i);
	  dVec ycoli;
	  _boundstrategy.to_f_representation(xcoli,ycoli);
	  ycandidates.col(i) = ycoli;
	}
      return ycandidates;
    }

    dVec pheno(const dVec &candidate)
    {
      // apply custom pheno function.
      dVec ncandidate;
      if (!_id)
	{
	  ncandidate = dVec(candidate.rows());
	  _phenof(candidate.data(),ncandidate.data(),candidate.rows());
	}

      // apply bounds.
      dVec phen = dVec::Zero(candidate.rows());
      if (_id)
	_boundstrategy.to_f_representation(candidate,phen);
      else _boundstrategy.to_f_representation(ncandidate,phen);

      // apply scaling.
      if (!_scalingstrategy._id)
	{
	  dVec sphen = dVec::Zero(phen.rows());
	  _scalingstrategy.scale_to_f(phen,sphen);
	  phen = sphen;
	}
      
      return phen;
    }
    
    dVec geno(const dVec &candidate)
    {
      dVec ccandidate = candidate;
      dVec gen = dVec::Zero(candidate.rows());
      
      // reverse scaling.
      if (!_scalingstrategy._id)
	{
	  _scalingstrategy.scale_to_internal(gen,candidate);
	  ccandidate = gen;
	}
      
      // reverse bounds.
      _boundstrategy.to_internal_representation(gen,ccandidate);

      // apply custom geno function.
      if (!_id)
	{
	  dVec ncandidate(gen.rows());
	  _genof(gen.data(),ncandidate.data(),gen.rows());
	  return ncandidate;
	}
      else return gen;
    }
    
    TBoundStrategy _boundstrategy;
    TransFunc _genof;
    TransFunc _phenof;
    bool _id; /**< geno/pheno transform is identity. */
    TScalingStrategy _scalingstrategy;
  };

  // specialization when no bound strategy nor scaling applies.
  template<> inline dMat GenoPheno<NoBoundStrategy,NoScalingStrategy>::pheno(const dMat &candidates)
    {
      if (_id)
	return candidates;
      else return pheno_candidates(candidates);
    }
  template<> inline dVec GenoPheno<NoBoundStrategy,NoScalingStrategy>::pheno(const dVec &candidate)
    {
      if (_id)
	return candidate;
      else
	{
	  dVec ncandidate(candidate.rows());
	  _phenof(candidate.data(),ncandidate.data(),candidate.rows());
	  return ncandidate;
	}
    }
  template<> inline dVec GenoPheno<NoBoundStrategy,NoScalingStrategy>::geno(const dVec &candidate)
    {
      if (_id)
	return candidate;
      else
	{
	  dVec ncandidate(candidate.rows());
	  _genof(candidate.data(),ncandidate.data(),candidate.rows());
	  return ncandidate;
	}
    }

  template<> inline dMat GenoPheno<NoBoundStrategy,linScalingStrategy>::pheno(const dMat &candidates)
    {
      dMat ncandidates;
      if (_id)
	ncandidates = candidates;
      else ncandidates = pheno_candidates(candidates);
      dMat ycandidates = dMat(ncandidates.rows(),ncandidates.cols());
#pragma omp parallel for if (candidates.cols() >= 100)
      for (int i=0;i<ncandidates.cols();i++)
	{
	  dVec xcoli = ncandidates.col(i);
	  dVec ycoli;
	  _scalingstrategy.scale_to_f(xcoli,ycoli);
	  ycandidates.col(i) = ycoli;
	}
      return ycandidates;
    }
}

#endif
