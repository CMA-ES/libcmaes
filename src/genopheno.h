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
#include <vector>

namespace libcmaes
{
  typedef std::function<void (const double*, double*, const int&)> TransFunc;
  
  template <class TBoundStrategy=NoBoundStrategy,class TScalingStrategy=NoScalingStrategy>
    class GenoPheno
    {
      friend class CMASolutions;
      
    public:
    GenoPheno()
    :_id(true)
    {}

    GenoPheno(TransFunc &genof, TransFunc &phenof)
    :_genof(genof),_phenof(phenof),_id(false)
    {}
    
    GenoPheno(const double *lbounds, const double *ubounds, const int &dim)
    :_boundstrategy(lbounds,ubounds,dim),_id(true),_scalingstrategy(lbounds,ubounds,dim)
    {
      if (_scalingstrategy._id)
	_boundstrategy = TBoundStrategy(lbounds,ubounds,dim);
      else
	{
	  std::vector<double> lb(dim,_scalingstrategy._intmin);
	  std::vector<double> ub(dim,_scalingstrategy._intmax);
	  _boundstrategy = TBoundStrategy(&lb.front(),&ub.front(),lbounds,ubounds,dim);
	}
    }

    GenoPheno(TransFunc &genof, TransFunc &phenof,
	      const double *lbounds, const double *ubounds, const int &dim)
    :_boundstrategy(lbounds,ubounds,dim),_genof(genof),_phenof(phenof),_id(false),_scalingstrategy(lbounds,ubounds,dim)
    {
      if (_scalingstrategy._id)
	_boundstrategy = TBoundStrategy(lbounds,ubounds,dim);
      else
	{
	  std::vector<double> lb(dim,_scalingstrategy._intmin);
	  std::vector<double> ub(dim,_scalingstrategy._intmax);
	  _boundstrategy = TBoundStrategy(&lb.front(),&ub.front(),lbounds,ubounds,dim);
	}
    }

    /**
     * \brief this is a dummy constructor to accomodate an easy to use 
     *        linear scaling with pwq bounds from a given scaling vector.
     *        Outside the library, the proper way to re-specialize for other
     *        custom scaling classes would be to inherit GenoPheno and 
     *        specialize constructors within the new class.
     * @param scaling vector for linear scaling of input parameters.
     */
    GenoPheno(const dVec &scaling,
	      const dVec &shift,
	      const double *lbounds=nullptr,
	      const double *ubounds=nullptr)
    :_id(true)
    {
      (void)scaling;
      (void)shift;
      (void)lbounds;
      (void)ubounds;
    }
    
    ~GenoPheno() {}

    private:
    dMat pheno_candidates(const dMat &candidates) const
    {
      if (!_id)
	{
	  dMat ncandidates = dMat(candidates.rows(),candidates.cols());
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

    dMat geno_candidates(const dMat &candidates) const
    {
      if (!_id)
	{
	  dMat ncandidates = dMat(candidates.rows(),candidates.cols());
#pragma omp parallel for if (candidates.cols() >= 100)
	  for (int i=0;i<candidates.cols();i++)
	    {
	      dVec in = dVec(candidates.rows());
	      _genof(candidates.col(i).data(),in.data(),candidates.rows());
	      ncandidates.col(i) = in;
	    }
	  return ncandidates;
	}
      return candidates;
    }
    
    public:
    dMat pheno(const dMat &candidates) const
    {
      // apply custom pheno function.
      dMat ncandidates = pheno_candidates(candidates);

      // apply bounds.
#pragma omp parallel for if (ncandidates.cols() >= 100)
      for (int i=0;i<ncandidates.cols();i++)
	{
	  dVec ycoli;
	  _boundstrategy.to_f_representation(ncandidates.col(i),ycoli);
	  ncandidates.col(i) = ycoli;
	}
      
      // apply scaling.
      if (!_scalingstrategy._id)
	{
#pragma omp parallel for if (ncandidates.cols() >= 100)
	  for (int i=0;i<ncandidates.cols();i++)
	    {
	      dVec ycoli;
	      _scalingstrategy.scale_to_f(ncandidates.col(i),ycoli);
	      ncandidates.col(i) = ycoli;
	    }
	}
      return ncandidates;
    }

    dMat geno(const dMat &candidates) const
    {
      // reverse scaling.
      dMat ncandidates = candidates;
      if (!_scalingstrategy._id)
	{
#pragma omp parallel for if (ncandidates.cols() >= 100)
	  for (int i=0;i<ncandidates.cols();i++)
	    {
	      dVec ycoli;
	      _scalingstrategy.scale_to_internal(ycoli,ncandidates.col(i));
	      ncandidates.col(i) = ycoli;
	    }
	}
      
      // reverse bounds.
#pragma omp parallel for if (ncandidates.cols() >= 100)
      for (int i=0;i<ncandidates.cols();i++)
	{
	  dVec ycoli;
	  _boundstrategy.to_internal_representation(ycoli,ncandidates.col(i));
	  ncandidates.col(i) = ycoli;
	}
      
      // apply custom geno function.
      ncandidates = geno_candidates(ncandidates);
      return ncandidates;
    }
    
    dVec pheno(const dVec &candidate) const
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
    
    dVec geno(const dVec &candidate) const
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

    TBoundStrategy get_boundstrategy() const { return _boundstrategy; }
      
    TBoundStrategy& get_boundstrategy_ref() { return _boundstrategy; }
      
    TScalingStrategy get_scalingstrategy() const { return _scalingstrategy; }

    void remove_dimensions(const std::vector<int> &k)
      {
	if (!_scalingstrategy.is_id())
	  _scalingstrategy.remove_dimensions(k);
     	if (!_boundstrategy.is_id())
	  _boundstrategy.remove_dimensions(k);
      }

    private:
    TBoundStrategy _boundstrategy;
    TransFunc _genof;
    TransFunc _phenof;
    bool _id = false; /**< geno/pheno transform is identity. */
    TScalingStrategy _scalingstrategy;
  };

  // specialization when no bound strategy nor scaling applies.
  template<> inline dMat GenoPheno<NoBoundStrategy,NoScalingStrategy>::pheno(const dMat &candidates) const
    {
      if (_id)
	return candidates;
      else return pheno_candidates(candidates);
    }
  template<> inline dVec GenoPheno<NoBoundStrategy,NoScalingStrategy>::pheno(const dVec &candidate) const
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
  template<> inline dVec GenoPheno<NoBoundStrategy,NoScalingStrategy>::geno(const dVec &candidate) const
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

  template<> inline dVec GenoPheno<NoBoundStrategy,linScalingStrategy>::pheno(const dVec &candidate) const
    {
      dVec ncandidate(candidate.rows());
      if (!_id)
	_phenof(candidate.data(),ncandidate.data(),candidate.rows());
      
      dVec sphen;
      if (!_id)
	_scalingstrategy.scale_to_f(ncandidate,sphen);
      else _scalingstrategy.scale_to_f(candidate,sphen);
      return sphen;
    }
  template<> inline dVec GenoPheno<NoBoundStrategy,linScalingStrategy>::geno(const dVec &candidate) const
    {
      dVec scand = dVec::Zero(candidate.rows());
      _scalingstrategy.scale_to_internal(scand,candidate);
      if (_id)
	return scand;
      else
	{
	  dVec ncandidate(candidate.rows());
	  _genof(scand.data(),scand.data(),candidate.rows());
	  return ncandidate;
	}
    }
  template<> inline dMat GenoPheno<NoBoundStrategy,linScalingStrategy>::pheno(const dMat &candidates) const
    {
      dMat ncandidates;
      if (!_id)
	ncandidates = pheno_candidates(candidates);
      else ncandidates = candidates;
      
      // apply scaling.
#pragma omp parallel for if (ncandidates.cols() >= 100)
      for (int i=0;i<ncandidates.cols();i++)
	{
	  dVec ycoli;
	  _scalingstrategy.scale_to_f(ncandidates.col(i),ycoli);
	  ncandidates.col(i) = ycoli;
	}
      return ncandidates;
    }
  
  template<> inline GenoPheno<NoBoundStrategy,linScalingStrategy>::GenoPheno(const dVec &scaling,
									     const dVec &shift,
									     const double *lbounds,
									     const double *ubounds)
    :_id(true)
    {
      (void)lbounds;
      (void)ubounds;
      _scalingstrategy = linScalingStrategy(scaling,shift);
    }
  
  template<> inline GenoPheno<pwqBoundStrategy,linScalingStrategy>::GenoPheno(const dVec &scaling,
									      const dVec &shift,
									      const double *lbounds,
									      const double *ubounds)
    :_id(true)
    {
      _scalingstrategy = linScalingStrategy(scaling,shift);
      if (lbounds == nullptr || ubounds == nullptr)
	return;
      dVec vlbounds = Eigen::Map<dVec>(const_cast<double*>(lbounds),scaling.size());
      dVec vubounds = Eigen::Map<dVec>(const_cast<double*>(ubounds),scaling.size());
      dVec nlbounds, nubounds;
      _scalingstrategy.scale_to_internal(nlbounds,vlbounds);
      _scalingstrategy.scale_to_internal(nubounds,vubounds);
      _boundstrategy = pwqBoundStrategy(nlbounds.data(),nubounds.data(),scaling.size());
    }
}

#endif
