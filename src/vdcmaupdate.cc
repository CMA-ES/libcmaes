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
 * Parts of this code are rewrittent from Y. Akimoto original code for
 * Y. Akimoto, A. Auger and N. Hansen, Comparison-Based Natural Gradient 
 * Optimization in High Dimension, GECCO-2014.
 */

#include "vdcmaupdate.h"
#include <iostream>

namespace libcmaes
{

  template <class TGenoPheno>
  void VDCMAUpdate::update(const CMAParameters<TGenoPheno> &parameters,
			   Eigen::EigenMultivariateNormal<double> &esolver,
			   CMASolutions &solutions)
  {
    (void)esolver; // esolver is not needed.
    
    // update of the mean
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * solutions._candidates.at(i).get_x_dvec();

    // reusable variables.
    dVec diffxmean = 1.0/solutions._sigma * (xmean-solutions._xmean); // (m^{t+1}-m^t)/sigma^t
    solutions._sepcsqinv = solutions._sepcov.cwiseInverse();
    
    // update psigma
    solutions._psigma = (1.0-parameters._csigma)*solutions._psigma;
    solutions._psigma += parameters._fact_ps * (dVec::Constant(parameters._dim,1.0) + solutions._v.cwiseProduct(solutions._v)).cwiseSqrt().cwiseInverse().cwiseProduct(solutions._sepcsqinv).cwiseProduct(diffxmean);
    double norm_ps = solutions._psigma.norm();
    
    // update pc
    solutions._hsig = 0;
    solutions._pc = (1.0-parameters._cc) * solutions._pc;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1-parameters._fixed_p.size()))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1;
    if (solutions._hsig)
      solutions._pc += parameters._fact_pc * diffxmean;
    
    // compute s and t
    double normv = solutions._v.squaredNorm(); //TODO: store in solution object.
    dVec vbar = solutions._v / std::sqrt(normv);
    dVec vbarbar = vbar.cwiseProduct(vbar);

    bool test = parameters._cmu + parameters._c1 * solutions._hsig > 0;
    if (test)
      {
	double gammav = 1.0+normv;
	double alpha = std::sqrt(normv*normv+(2-1.0/std::sqrt(gammav))*gammav / (vbarbar.maxCoeff())) / (2.0+normv); // Eq. (7)
	alpha = std::min(1.0,alpha);
	double b = -(1-alpha*alpha)*normv*normv/gammav + 2.0*alpha*alpha; // Lemma 3.4
	dVec Ainv = (dVec::Constant(parameters._dim,2.0) - (b + 2*alpha*alpha) * vbarbar).cwiseInverse(); // paper has alpha, code has alpha^2
	dVec Ainvbb = Ainv.cwiseProduct(vbarbar);
	
	dMat y(parameters._dim,parameters._mu);
	for (int i=0;i<parameters._mu;i++)
	  {
	    y.col(i) = solutions._sepcsqinv.cwiseProduct(solutions._candidates.at(i).get_x_dvec() - solutions._xmean) / solutions._sigma;
	  }
        dMat ym = solutions._sepcsqinv.cwiseProduct(solutions._pc);
	dVec yvbar = vbar.transpose() * y;
	dVec pvec = (y.cwiseProduct(y) - (normv / gammav) * ((vbar*yvbar.transpose()).cwiseProduct(y)) - dMat::Constant(parameters._dim,parameters._mu,1.0))*parameters._weights;
	dMat qvec = (y.cwiseProduct(yvbar.replicate(1,parameters._dim).transpose()) - 0.5*vbar*(yvbar.cwiseProduct(yvbar) + dVec::Constant(parameters._mu,gammav)).transpose())*parameters._weights;
	yvbar = vbar.transpose() * ym;
	dVec pone = (ym.cwiseProduct(ym) - (normv / gammav) * ((vbar*yvbar.transpose()).cwiseProduct(ym)) - dVec::Constant(parameters._dim,1.0));
	dVec qone = (ym.cwiseProduct(yvbar.replicate(1,parameters._dim).transpose()) - 0.5*vbar*(yvbar.cwiseProduct(yvbar) + dVec::Constant(1,gammav)));
	
	if (solutions._hsig)
	  {
	    pvec = parameters._cmu * pvec + parameters._c1 * pone;
	    qvec = parameters._cmu * qvec + parameters._c1 * qone;
	  }
	else
	  {
	    pvec = parameters._cmu * pvec;
	    qvec = parameters._cmu * qvec;
	  }
	
	double nu = (vbar.transpose()*qvec)(0);
	dVec rvec = pvec - (alpha/gammav) * ((2.0+normv)*(vbar.cwiseProduct(qvec)) - normv*nu*vbarbar);
	nu = Ainvbb.transpose()*rvec;
	double nu2 = Ainvbb.transpose()*vbarbar;
	dVec svec = rvec.cwiseProduct(Ainv) - (b * nu/(1.0+b*nu2)) * Ainvbb;
	nu = svec.transpose()*vbarbar;
	dVec ngv = (qvec - alpha * ((2+normv) * (vbar.cwiseProduct(svec)) - nu * vbar))/std::sqrt(normv);
	dVec ngd = solutions._sepcov.cwiseProduct(svec);
	
	solutions._v += ngv;
	solutions._sepcov += ngd;
      }

    // update sigma.
    if (parameters._tpa < 2)
      solutions._sigma *= exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0));
    else if (solutions._niter > 0)
      solutions._sigma *= std::exp(solutions._tpa_s / parameters._dsigma);
    
    // set mean.
    if (parameters._tpa)
      solutions._xmean_prev = solutions._xmean;
    solutions._xmean = xmean;
  }

  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  template void VDCMAUpdate::update(const CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>&,Eigen::EigenMultivariateNormal<double>&,CMASolutions&);
  
}
