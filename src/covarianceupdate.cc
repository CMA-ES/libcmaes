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

#include "covarianceupdate.h"

namespace libcmaes
{

  void CovarianceUpdate::update(const CMAParameters &parameters,
				EigenMultivariateNormal<double> &esolver,
				CMASolutions &solutions)
  {
    // compute mean, Eq. (2)
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * solutions._candidates.at(i)._x;
    
    // reusable variables.
    dVec diffxmean = 1.0/solutions._sigma * (xmean-solutions._xmean); // (m^{t+1}-m^t)/sigma^t
    if (solutions._updated_eigen)
      solutions._csqinv = esolver._eigenSolver.operatorInverseSqrt();

    // update psigma, Eq. (3)
    solutions._psigma = (1.0-parameters._csigma)*solutions._psigma
      + parameters._fact_ps * solutions._csqinv * diffxmean;
    double norm_ps = solutions._psigma.norm();

    // update pc, Eq. (4)
    solutions._hsig = 0;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1; //TODO: simplify equation instead.
    solutions._pc = (1.0-parameters._cc) * solutions._pc + solutions._hsig * parameters._fact_pc * diffxmean;
    dMat spc = solutions._pc * solutions._pc.transpose();
    
    // covariance update, Eq (5).
    dMat wdiff = dMat::Zero(parameters._dim,parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec difftmp = solutions._candidates.at(i)._x - solutions._xmean;
	wdiff += parameters._weights[i] * (difftmp*difftmp.transpose());
      }
    wdiff *= 1.0/(solutions._sigma*solutions._sigma);
    solutions._cov = (1-parameters._c1-parameters._cmu+(1-solutions._hsig)*parameters._c1*parameters._cc*(2.0-parameters._cc))*solutions._cov + parameters._c1*spc + parameters._cmu*wdiff;
    
    // sigma update, Eq. (6)
    solutions._sigma *= exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0));

    // set mean.
    solutions._xmean = xmean;
  }
  
}
