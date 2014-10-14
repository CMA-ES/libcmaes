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

#include "vdcmaupdate.h"
#include <iostream>

namespace libcmaes
{

  template <class TGenoPheno>
  void VDCMAUpdate::update(const CMAParameters<TGenoPheno> &parameters,
			   EigenMultivariateNormal<double> &esolver,
			   CMASolutions &solutions)
  {
    (void)esolver; // esolver is not needed.
    
    // update of the mean
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * solutions._candidates.at(i).get_x_dvec(); //TODO: beware

    // reusable variables.
    dVec diffxmean = 1.0/solutions._sigma * (xmean-solutions._xmean); // (m^{t+1}-m^t)/sigma^t
    solutions._sepcsqinv = solutions._sepcov.cwiseInverse();
    
    // update psigma
    solutions._psigma = (1.0-parameters._csigma)*solutions._psigma;
    solutions._psigma += parameters._fact_ps * ((dVec::Constant(parameters._dim,1.0) + solutions._v.cwiseProduct(solutions._v)).cwiseSqrt().cwiseInverse().cwiseProduct(solutions._sepcsqinv)).cwiseProduct(diffxmean);
    double norm_ps = solutions._psigma.norm();
    
    // update pc
    solutions._hsig = 0;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1-parameters._fixed_p.size()))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1; //TODO: simplify equation instead.
    solutions._pc = (1.0-parameters._cc) * solutions._pc + solutions._hsig * parameters._fact_pc * diffxmean;
    dMat spc = solutions._pc.cwiseProduct(solutions._pc);

    // sigma update
    //solutions._sigma *= exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0));
    
    // compute s and t
    double normv = solutions._v.squaredNorm(); //TODO: store in solution object.
    dVec vbar = solutions._v / std::sqrt(normv);
    dVec vbarbar = vbar.cwiseProduct(vbar);
    double gammav = 1.0+normv;
    double alpha = std::sqrt(normv*normv+(2-1.0/std::sqrt(gammav))*gammav / (vbarbar.maxCoeff())) / (2.0+normv); // Eq. (7)
    alpha = std::min(1.0,alpha);
    double b = -(1-alpha*alpha)*normv*normv/gammav + 2.0*alpha*alpha; // Lemma 3.4
    dVec Ainv = (dVec::Constant(parameters._dim,2.0) - (b + 2*alpha) * vbarbar).cwiseInverse();

    dVec nsepcov = solutions._sepcov;
    dVec s(parameters._dim);
    dVec t(parameters._dim);
    //TODO: parallelize
    for (int i=0;i<parameters._mu+1;i++)
      {
	dVec y;
	if (i != parameters._mu)
	  y = solutions._sepcsqinv.cwiseProduct(solutions._candidates.at(i).get_x_dvec() - solutions._xmean) / solutions._sigma;
	else y = solutions._sepcsqinv.cwiseProduct(solutions._pc);
	double yvbar = y.dot(vbar);
	s = y.cwiseProduct(y) - (normv*yvbar/gammav) * y.cwiseProduct(vbar) - dVec::Constant(parameters._dim,1.0); // step 1
	t = yvbar*y - 0.5*(yvbar*yvbar + gammav)*vbar;
	s -= alpha/gammav * ((2+normv)*vbar.cwiseProduct(t) - normv*vbar.dot(t)*vbarbar);
	s = Ainv.cwiseProduct(s) - (b/(1+b*vbarbar.dot(Ainv.cwiseProduct(vbarbar)))) * s.dot(Ainv.cwiseProduct(vbarbar))*Ainv.cwiseProduct(vbarbar);
	t -= alpha*((2+normv)*vbar.cwiseProduct(s)-s.dot(vbarbar)*vbar);
	
	// covariance update.
	if (i != parameters._mu)
	  {
	    solutions._v += parameters._cmu * (parameters._weights[i]/std::sqrt(normv)) * t;
	    nsepcov += parameters._cmu * parameters._weights[i] * solutions._sepcov.cwiseProduct(s);
	  }
	else
	  {
	    solutions._v += (1-solutions._hsig)*parameters._c1*t/std::sqrt(normv);
	    nsepcov += (1-solutions._hsig)*parameters._c1*solutions._sepcov.cwiseProduct(s);
	  }
      }
    solutions._sepcov = nsepcov;

    solutions._sigma *= exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0));
    
    // set mean.
    solutions._xmean = xmean;
  }

  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy>>&,EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy>>&,EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>&,EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>&,EigenMultivariateNormal<double>&,CMASolutions&);
  
}
