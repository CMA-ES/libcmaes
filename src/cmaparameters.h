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

#ifndef CMAPARAMETERS_H
#define CMAPARAMETERS_H

#include "parameters.h"
#include "eo_matrix.h"
#include <float.h>

namespace libcmaes
{
  class CMAParameters : public Parameters
  {
  public:
    CMAParameters() {}; //TODO: var init even if this constructor is not supposed to be used for now.
    CMAParameters(const int &dim, const int &lambda=-1,
		  const int &max_iter=-1, const int &max_fevals=-1,
		  const std::string &fplot="",
		  const double &sigma_init=-1.0,
		  const double &x0=std::numeric_limits<double>::min(),
		  const uint64_t &seed=0);
    ~CMAParameters();
    
    int _mu;
    dVec _weights;
    double _csigma;
    double _c1;
    double _cmu;
    double _cc;
    double _muw;
    double _dsigma;
    
    // computed once at init for speeding up operations.
    double _fact_ps;
    double _fact_pc;
    double _chi; // norm of N(0,I).

    double _sigma_init;

    int _nrestarts; // when applicable.
    bool _lazy_update; /**< covariance lazy update. */
    double _lazy_value; /**< reference trigger for lazy update. */

    // active cma.
    double _cm; /**< learning rate for the mean. */
    double _alphacov;
    double _alphaminusold;
    double _deltamaxsigma;
    double _lambdamintarget;
    double _alphaminusmin;
  };
  
}

#endif
