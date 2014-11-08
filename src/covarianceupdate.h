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

#ifndef COVARIANCEUPDATE_H
#define COVARIANCEUPDATE_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include "eigenmvn.h"

namespace libcmaes
{

  /**
   * \brief Covariance Matrix update.
   *        This is an implementation closely follows:
   * Hansen, N. (2009). Benchmarking a BI-Population CMA-ES on the BBOB-2009 Function Testbed. Workshop Proceedings of the GECCO Genetic and Evolutionary Computation Conference, ACM, pp. 2389-2395
   */
  class CovarianceUpdate
  {
  public:
    /**
     * \brief update the covariance matrix.
     * @param parameters current set of parameters
     * @param esolver Eigen eigenvalue solver
     * @param solutions currrent set of solutions.
     */
    template <class TGenoPheno>
    static void update(const CMAParameters<TGenoPheno> &parameters,
		       Eigen::EigenMultivariateNormal<double> &esolver,
		       CMASolutions &solutions);
  };
  
}

#endif
