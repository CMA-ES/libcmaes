
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
