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

#include "acovarianceupdate.h"
#include <iostream>

namespace libcmaes
{

  template <class TGenoPheno>
  void ACovarianceUpdate::update(const CMAParameters<TGenoPheno> &parameters,
				 Eigen::EigenMultivariateNormal<double> &esolver,
				 CMASolutions &solutions)
  {
    // compute mean, Eq. (2)
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * (solutions._candidates.at(i).get_x_dvec() - solutions._xmean);
    xmean *= parameters._cm;
    xmean += solutions._xmean;
  
     // reusable variables.
    dVec diffxmean = 1.0/(solutions._sigma*parameters._cm) * (xmean-solutions._xmean); // (m^{t+1}-m^t)/(c_m*sigma^t)
    if (solutions._updated_eigen && !parameters._sep)
      solutions._csqinv = esolver._eigenSolver.operatorInverseSqrt();
    else if (parameters._sep)
      solutions._sepcsqinv = solutions._sepcov.cwiseInverse().cwiseSqrt();
    
    // update psigma, Eq. (3)
    if (!parameters._sep)
      solutions._psigma = (1.0-parameters._csigma)*solutions._psigma
	+ parameters._fact_ps * solutions._csqinv * diffxmean;
    else solutions._psigma = (1.0-parameters._csigma)*solutions._psigma
	   + parameters._fact_ps * solutions._sepcsqinv.cwiseProduct(diffxmean);
    double norm_ps = solutions._psigma.norm();

    // update pc, Eq. (4-5)
    solutions._hsig = 0;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1-parameters._fixed_p.size()))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1; //TODO: simplify equation instead.
    solutions._pc = (1.0-parameters._cc) * solutions._pc + solutions._hsig * parameters._fact_pc * diffxmean;
    dMat spc;
    if (!parameters._sep)
      spc = solutions._pc * solutions._pc.transpose();
    else spc = solutions._pc.cwiseProduct(solutions._pc);
    
    // Cmu+, Eq. (6)
    dMat cmuplus;
    if (!parameters._sep)
      cmuplus = dMat::Zero(parameters._dim,parameters._dim);
    else cmuplus = dMat::Zero(parameters._dim,1);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec difftmp = solutions._candidates.at(i).get_x_dvec() - solutions._xmean;
	if (!parameters._sep)
	  cmuplus += parameters._weights[i] * (difftmp*difftmp.transpose());
	else cmuplus += parameters._weights[i] * (difftmp.cwiseProduct(difftmp));
      }
    cmuplus *= 1.0/(solutions._sigma*solutions._sigma);
        
    // Cmu-, Eq. (7)
    dMat cmuminus;
    if (!parameters._sep)
      cmuminus = dMat::Zero(parameters._dim,parameters._dim);
    else cmuminus = dMat::Zero(parameters._dim,1);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec ytmp = solutions._candidates.at(parameters._lambda-i-1).get_x_dvec()-solutions._xmean;
	//dVec yl = (solutions._csqinv * (solutions._candidates.at(parameters._lambda-parameters._mu+i)._x-solutions._xmean)).norm() / (solutions._csqinv * ytmp).norm() * ytmp * 1.0/solutions._sigma;
	dVec yl = ytmp * 1.0/solutions._sigma; // NH says this is a good enough value.
	if (!parameters._sep)
	  cmuminus += parameters._weights[i] * yl*yl.transpose();
	else cmuminus += parameters._weights[i] * yl.cwiseProduct(yl);
      }
    
    // covariance update, Eq. (8)
    dMat cminusdenom;
    if (!parameters._sep)
      cminusdenom = solutions._csqinv*cmuminus*solutions._csqinv;
    else cminusdenom = solutions._sepcsqinv.cwiseProduct(cmuminus.cwiseProduct(solutions._sepcsqinv));
    double cminustmp = parameters._lambdamintarget;
    if (!parameters._sep)
      {
	Eigen::SelfAdjointEigenSolver<dMat> tmpesolve(cminusdenom); // XXX: computing eigenvalues, could be avoid with an upper bound.
	cminustmp = tmpesolve.eigenvalues().maxCoeff();
      }
    else cminustmp = cminusdenom.maxCoeff();
    double cminusmin = parameters._alphaminusmin * (1.0-parameters._cmu)*(1.0-parameters._lambdamintarget) / cminustmp;
    double cminus = std::min(cminusmin,(1-parameters._cmu)*parameters._alphacov/8.0*(parameters._muw/(pow(parameters._dim+2.0,1.5)+2.0*parameters._muw)));
    if (!parameters._sep)
      solutions._cov = (1-parameters._c1-parameters._cmu + cminus*parameters._alphaminusold)*solutions._cov + parameters._c1*spc + (parameters._cmu + cminus * (1.0-parameters._alphaminusold))*cmuplus - cminus * cmuminus;
    else solutions._sepcov = (1-parameters._c1-parameters._cmu + cminus*parameters._alphaminusold)*solutions._sepcov + parameters._c1*spc + (parameters._cmu + cminus * (1.0-parameters._alphaminusold))*cmuplus - cminus * cmuminus;
    
    // sigma update, Eq. (9)
    if (parameters._tpa < 2)
      solutions._sigma *= std::min(parameters._deltamaxsigma,exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0)));
    else if (solutions._niter > 0)
      solutions._sigma *= std::exp(solutions._tpa_s / parameters._dsigma);
        
    // set mean.
    solutions._xmean_prev = solutions._xmean;
    solutions._xmean = xmean;
  }

  template void ACovarianceUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void ACovarianceUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void ACovarianceUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void ACovarianceUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
}
