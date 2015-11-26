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
#include <vector>

namespace libcmaes
{

  /**
   * \brief profile likelihood object holder as a set of points and values.
   */
  class pli
  {
    friend class CMASolutions;
    template <class U> friend class errstats;
    
  public:
    pli() {}
    
    /**
     * \brief profile likelihood constructor
     * @param k dimension in which the profile likelihood was computed
     * @param samplesize number of steps of the linesearch direction
     * @param dim dimension of the objective function parameter space
     * @param xm vector of parameters at fvalue
     * @param fvalue the function minima around which the profile likelihood was computed
     * @param fup the function deviation for which the profile likelihood was computed
     * @param delta tolerance around fvalue + fup for which the profile likelihood was computed
     */
  pli(const int &k, const int &samplesize, const int &dim,
      const dVec &xm, const double &fvalue, const double &fup, const double &delta)
    :_k(k),_samplesize(samplesize),_fvaluem(dVec::Zero(2*samplesize+1)),_xm(dMat::Zero(2*samplesize+1,dim)),_min(0.0),_max(0.0),_err(2*samplesize+1),_fup(fup),_delta(delta)
      {
	_fvaluem[samplesize] = fvalue;
	_xm.row(samplesize) = xm.transpose();
	_err[samplesize] = 1; // should be current sol status...
      }
    
    ~pli() {};
    
    /**
     * \brief find bounds around the objective function parameters for a given value of f,
     *        base on pre-computed profile likelihood points.
     * @param fvalue function value
     * @param minindex index of the profile likelihood point that is the lower bound
     * @param maxindex index of the profile likelihood point that is the upper bound
     */
    std::pair<double,double> getMinMax(const double &fvalue,
				       int &minindex, int &maxindex)
    {
      (_fvaluem.head(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&minindex);
      (_fvaluem.tail(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&maxindex);
      double min = _xm(minindex,_k);
      double max = _xm(_samplesize + 1 + maxindex,_k);
      if (min > max)
	std::swap(min,max);
      return std::pair<double,double>(min,max);
    }

    /**
     * \brie sets the bounds for this profile likelihood object based on original function value + fup
     */
    void setMinMax()
    {
      std::pair<double,double> mm = getMinMax(_fvaluem[_samplesize]+_fup,_minindex,_maxindex);
      _min = mm.first;
      _max = mm.second;
    }

    /**
     * \brief sets the errors bounds for this profile likelihood.
     */
    void setErrMinMax()
    {
      setMinMax();
      _errmin = _min - _xm(_samplesize,_k);
      _errmax = _max - _xm(_samplesize,_k);
    }

    /**
     * \brief get lower error bound
     * @return lower error bound
     */
    inline double get_err_min() const
    {
      return _errmin;
    }

    /**
     * \brief get upper error bound
     * @return upper error bound
     */
    inline double get_err_max() const
    {
      return _errmax;
    }
    
    // accessors
    inline int get_k() const
    {
      return _k;
    }

    inline int get_samplesize() const
    {
      return _samplesize;
    }
    
    inline dVec get_fvaluem() const
    {
      return _fvaluem;
    }

    inline dMat get_xm() const
    {
      return _xm;
    }

    inline double get_min() const
    {
      return this->_min;
    }

    inline double get_max() const
    {
      return _max;
    }

    private:
    int _k = -1;
    int _samplesize = 0;
    dVec _fvaluem;
    dMat _xm; // in phenotype
    double _min = 0.0;
    double _max = 0.0;
    double _errmin = 0.0;
    double _errmax = 0.0;
    int _minindex = -1;
    int _maxindex = -1;
    std::vector<int> _err; // errors from profile likelihood computations as run status codes.
    double _fup; // the function deviation for which this profile likelihood was computed.
    double _delta; // the tolerance around fvalue + fup for which this profile likelihood was computed.
  };

}

#endif
