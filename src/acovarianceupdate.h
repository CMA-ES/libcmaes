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

#ifndef ACOVARIANCEUPDATE_H
#define ACOVARIANCEUPDATE_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include "eigenmvn.h"

namespace libcmaes
{

  /**
   * \brief Active Covariance Matrix update.
   *        This implementation closely follows
   *        N. Hansen, R. Ros, "Benchmarking a Weighted Negative Covariance Matrix 
   *                            Update on the BBOB-2010 Noiseless Testbed", GECCO'10, 2010.
   */
  class CMAES_EXPORT ACovarianceUpdate
  {
  public:
    template <class TGenoPheno>
    static void update(const CMAParameters<TGenoPheno> &parameters,
		       Eigen::EigenMultivariateNormal<double> &esolver,
		       CMASolutions &solutions);
  };
  
}

#endif
