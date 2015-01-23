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

#include "covarianceupdate.h"
#include <iostream>

namespace libcmaes
{

  template <class TGenoPheno>
  void CovarianceUpdate::update(const CMAParameters<TGenoPheno> &parameters,
				Eigen::EigenMultivariateNormal<double> &esolver,
				CMASolutions &solutions)
  {
    // compute mean, Eq. (2)
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * solutions._candidates.at(i).get_x_dvec();
    
    // reusable variables.
    dVec diffxmean = 1.0/solutions._sigma * (xmean-solutions._xmean); // (m^{t+1}-m^t)/sigma^t
    if (solutions._updated_eigen && !parameters._sep) //TODO: shall not recompute when using gradient, as it is computed in ask.
      solutions._csqinv = esolver._eigenSolver.operatorInverseSqrt();
    else if (parameters._sep)
      solutions._sepcsqinv = solutions._sepcov.cwiseInverse().cwiseSqrt();
    
    // update psigma, Eq. (3)
    solutions._psigma = (1.0-parameters._csigma)*solutions._psigma;
    if (!parameters._sep)
      solutions._psigma += parameters._fact_ps * solutions._csqinv * diffxmean;
    else
      solutions._psigma += parameters._fact_ps * solutions._sepcsqinv.cwiseProduct(diffxmean);
    double norm_ps = solutions._psigma.norm();

    // update pc, Eq. (4)
    solutions._hsig = 0;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1-parameters._fixed_p.size()))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1; //TODO: simplify equation instead.
    solutions._pc = (1.0-parameters._cc) * solutions._pc + solutions._hsig * parameters._fact_pc * diffxmean;
    dMat spc;
    if (!parameters._sep)
      spc = solutions._pc * solutions._pc.transpose();
    else spc = solutions._pc.cwiseProduct(solutions._pc);
    
    // covariance update, Eq (5).
    dMat wdiff;
    if (!parameters._sep)
      wdiff = dMat::Zero(parameters._dim,parameters._dim);
    else wdiff = dMat::Zero(parameters._dim,1);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec difftmp = solutions._candidates.at(i).get_x_dvec() - solutions._xmean;
	if (!parameters._sep)
	  wdiff += parameters._weights[i] * (difftmp*difftmp.transpose());
	else wdiff += parameters._weights[i] * (difftmp.cwiseProduct(difftmp));
      }
    wdiff *= 1.0/(solutions._sigma*solutions._sigma);
    if (!parameters._sep)
      solutions._cov = (1-parameters._c1-parameters._cmu+(1-solutions._hsig)*parameters._c1*parameters._cc*(2.0-parameters._cc))*solutions._cov + parameters._c1*spc + parameters._cmu*wdiff;
    else
      {
	solutions._sepcov = (1-parameters._c1-parameters._cmu+(1-solutions._hsig)*parameters._c1*parameters._cc*(2.0-parameters._cc))*solutions._sepcov + parameters._c1*spc + parameters._cmu*wdiff;
      }
    
    // sigma update, Eq. (6)
    if (parameters._tpa < 2)
      solutions._sigma *= std::exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0));
    else if (solutions._niter > 0)
      solutions._sigma *= std::exp(solutions._tpa_s / parameters._dsigma);
    
    // set mean.
    if (parameters._tpa)
      solutions._xmean_prev = solutions._xmean;
    solutions._xmean = xmean;
  }

  template void CovarianceUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void CovarianceUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void CovarianceUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void CovarianceUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
}
