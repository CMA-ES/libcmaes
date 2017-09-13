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

#ifndef VDCMAUPDATE_H
#define VDCMAUPDATE_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include "eigenmvn.h"

namespace libcmaes
{
  /**
   * \brief VD-CMA update that is a linear time/space variant of CMA-ES
   * This is an implementation that closely follows:
   * Y. Akimoto, A. Auger and N. Hansen: Comparison-Based Natural 
   *  Gradient Optimization in High Dimension. In Proceedings of Genetic
   *  and Evolutionary Computation Conference (2014)
   */
  class CMAES_EXPORT VDCMAUpdate
  {
  public:
    template <class TGenoPheno>
      static void update(const CMAParameters<TGenoPheno> &parameters,
			 Eigen::EigenMultivariateNormal<double> &esolver,
			 CMASolutions &solutions);
  };
  
}

#endif
