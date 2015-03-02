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

/**
 * This is based on Nikolaus Hansen's code from the C version of CMA-ES,
 * see boundary_transformation.c from the package available from
 * https://www.lri.fr/~hansen/cmaes_inmatlab.html#C
 */

#include "pwq_bound_strategy.h"
#include "llogging.h"
#include <limits>
#include <iostream>

namespace libcmaes
{
  pwqBoundStrategy::pwqBoundStrategy()
  {
  }
  
  pwqBoundStrategy::pwqBoundStrategy(const double *lbounds, const double *ubounds, const int &dim)
    :_lbounds(Eigen::Map<dVec>(const_cast<double*>(lbounds),dim)),_ubounds(Eigen::Map<dVec>(const_cast<double*>(ubounds),dim)),
     _phenolbounds(_lbounds),_phenoubounds(_ubounds)
  {
    // init al & ul.
    dVec tmpdiff1 = _ubounds - _lbounds;
    dVec tmpdiff2 = 0.5*tmpdiff1;
    dVec tmpal = (1.0/20.0) * (dVec::Constant(dim,1.0) + _lbounds.cwiseAbs());
    _al = tmpdiff2.cwiseMin(tmpal);
    
    dVec tmpau = (1.0/20.0) * (dVec::Constant(dim,1.0) + _ubounds.cwiseAbs());
    _au = tmpdiff2.cwiseMin(tmpau);
    
    // compute static variables.
    _xlow = _lbounds - 2.0 * _al - tmpdiff2;
    _xup = _ubounds + 2.0 * _au + tmpdiff2;
    _r = 2.0 * (tmpdiff1 + _al + _au);
  }

  pwqBoundStrategy::pwqBoundStrategy(const double *lbounds, const double *ubounds,
				     const double *plbounds, const double *pubounds, const int &dim)
    :_lbounds(Eigen::Map<dVec>(const_cast<double*>(lbounds),dim)),_ubounds(Eigen::Map<dVec>(const_cast<double*>(ubounds),dim)),
     _phenolbounds(Eigen::Map<dVec>(const_cast<double*>(plbounds),dim)),_phenoubounds(Eigen::Map<dVec>(const_cast<double*>(pubounds),dim))
  {
    // init al & ul.
    dVec tmpdiff1 = _ubounds - _lbounds;
    dVec tmpdiff2 = 0.5*tmpdiff1;
    dVec tmpal = (1.0/20.0) * (dVec::Constant(dim,1.0) + _lbounds.cwiseAbs());
    _al = tmpdiff2.cwiseMin(tmpal);
    
    dVec tmpau = (1.0/20.0) * (dVec::Constant(dim,1.0) + _ubounds.cwiseAbs());
    _au = tmpdiff2.cwiseMin(tmpau);
    
    // compute static variables.
    _xlow = _lbounds - 2.0 * _al - tmpdiff2;
    _xup = _ubounds + 2.0 * _au + tmpdiff2;
    _r = 2.0 * (tmpdiff1 + _al + _au);
  }
  
  pwqBoundStrategy::~pwqBoundStrategy()
  {
  }

  void pwqBoundStrategy::to_f_representation(const dVec &x, dVec &y) const
  {
    shift_into_feasible(x,y);
    for (int i=0;i<x.rows();i++) //TODO: vectorize ?
      {
	if (y[i] < _lbounds[i] + _al[i])
	  y[i] = _lbounds[i] + (y[i] - (_lbounds[i] - _al[i])) * (y[i] - (_lbounds[i] - _al[i])) / 4.0 / _al[i];
	else if (y[i] > _ubounds[i] - _au[i])
	  y[i] = _ubounds[i] - (y[i] - (_ubounds[i] + _au[i])) * (y[i] - (_ubounds[i] + _au[i])) / 4.0 / _au[i];
      }
  }

  void pwqBoundStrategy::shift_into_feasible(const dVec &x, dVec &x_s) const
  {
    x_s = x;
    for (int i=0;i<x.rows();i++) //TODO: vectorize ?
      {
	if (x_s[i] < _xlow[i])
	  x_s[i] += _r[i] * (1.0 + static_cast<int>((_xlow[i]-x_s[i])/_r[i])); // shift up.
	if (x_s[i] > _xup[i])
	  x_s[i] -= _r[i] * (1.0 + static_cast<int>((x_s[i]-_xup[i])/_r[i])); // shift down;
	if (x_s[i] < _lbounds[i] - _al[i])
	  x_s[i] += 2.0 * (_lbounds[i] - _al[i] - x_s[i]);
	if (x_s[i] > _ubounds[i] + _au[i])
	  x_s[i] -= 2.0 * (x_s[i] - _ubounds[i] - _au[i]);

	if ((x_s[i] < _lbounds[i] - _al[i] - 1e-15) || (x_s[i] > _ubounds[i] + _au[i] + 1e-15))
	  {
	    LOG(FATAL) << "error in shifting pwq bounds in dimension " << i << ": lb=" << _lbounds[i] << " / ub=" << _ubounds[i] << " / al=" << _al[i] << " / au=" << _au[i] << " / x_s=" << x_s[i] << " / x=" << x[i] << " / xlow=" << _xlow[i] << " / xup=" << _xup[i] << " / r=" << _r[i] << std::endl;
	  }
      }
  }

  void pwqBoundStrategy::to_internal_representation(dVec &x, const dVec &y) const
  {
    x = y;
    for (int i=0;i<y.rows();i++)
      {
	if (x[i] < _lbounds[i] + _al[i])
	  x[i] = (_lbounds[i] - _al[i]) + 2.0 * sqrt(_al[i] * fabs(_lbounds[i] - x[i]));
	else if (x[i] > _ubounds[i] - _au[i])
	  x[i] = (_ubounds[i] + _au[i]) - 2.0 * sqrt(_au[i] * fabs(x[i] - _ubounds[i]));
      }
  }

  void pwqBoundStrategy::remove_dimensions(const std::vector<int> &k)
  {
    for (const int i: k)
      {
	removeElement(_lbounds,i);
	removeElement(_ubounds,i);
	removeElement(_al,i);
	removeElement(_au,i);
	removeElement(_xlow,i);
	removeElement(_xup,i);
	removeElement(_r,i);
	removeElement(_phenolbounds,i);
	removeElement(_phenoubounds,i);
      }
  }
  
}
