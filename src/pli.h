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

#ifndef PLI_H
#define PLI_H

#include "eo_matrix.h"

namespace libcmaes
{

  /**
   * \brief profile likelihood as a set of points and values. 
   */
  class pli
  {
  public:
  pli() {}
  pli(const int &k, const int &samplesize, const int &dim,
      const dVec &xm, const double &fvalue)
    :_k(k),_samplesize(samplesize),_fvaluem(dVec::Zero(2*samplesize+1)),_xm(dMat::Zero(2*samplesize+1,dim)),_min(0.0),_max(0.0)
      {
	_fvaluem[samplesize] = fvalue;
	_xm.row(samplesize) = xm.transpose();
      }
    ~pli() {};

    void setMinMax()
    {
      _min = _xm(0,_k);
      _max = _xm(2*_samplesize,_k);
      if (_min > _max)
	std::swap(_min,_max);
    }

    void setErrMinMax()
    {
      setMinMax();
      _errmin = _min - _xm(_samplesize,_k);
      _errmax = _max - _xm(_samplesize,_k);
    }
    
    std::pair<double,double> getMinMax(const double &fvalue)
    {
      dMat::Index mindex[2];
      (_fvaluem.head(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&mindex[0]);
      (_fvaluem.tail(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&mindex[1]);
      double min = _xm(mindex[0],_k);
      double max = _xm(_samplesize + 1 + mindex[1],_k);
      if (min > max)
	std::swap(min,max);
      return std::pair<double,double>(min,max);
    }
    
    int _k = -1;
    int _samplesize = 0;
    dVec _fvaluem;
    dMat _xm;
    double _min = 0.0;
    double _max = 0.0;
    double _errmin = 0.0;
    double _errmax = 0.0;
  };

}

#endif
