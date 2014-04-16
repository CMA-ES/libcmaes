
#include "acovarianceupdate.h"
#include <iostream>

namespace libcmaes
{

  void ACovarianceUpdate::update(const CMAParameters &parameters,
				 EigenMultivariateNormal<double> &esolver,
				 CMASolutions &solutions)
  {
    // compute mean, Eq. (2)
    dVec xmean = dVec::Zero(parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      xmean += parameters._weights[i] * (solutions._candidates.at(i)._x - solutions._xmean);
    xmean *= parameters._cm;
    xmean += solutions._xmean;
  
     // reusable variables.
    dVec diffxmean = 1.0/(solutions._sigma*parameters._cm) * (xmean-solutions._xmean); // (m^{t+1}-m^t)/(c_m*sigma^t)
    if (solutions._updated_eigen)
      solutions._csqinv = esolver._eigenSolver.operatorInverseSqrt();

    // update psigma, Eq. (3)
    solutions._psigma = (1.0-parameters._csigma)*solutions._psigma
      + parameters._fact_ps * solutions._csqinv * diffxmean;
    double norm_ps = solutions._psigma.norm();

    // update pc, Eq. (4-5)
    solutions._hsig = 0;
    double val_for_hsig = sqrt(1.0-pow(1.0-parameters._csigma,2.0*(solutions._niter+1)))*(1.4+2.0/(parameters._dim+1))*parameters._chi;
    if (norm_ps < val_for_hsig)
      solutions._hsig = 1; //TODO: simplify equation instead.
    solutions._pc = (1.0-parameters._cc) * solutions._pc + solutions._hsig * parameters._fact_pc * diffxmean;
    dMat spc = solutions._pc * solutions._pc.transpose();

    // Cmu+, Eq. (6)
    dMat cmuplus = dMat::Zero(parameters._dim,parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec difftmp = solutions._candidates.at(i)._x - solutions._xmean;
	cmuplus += parameters._weights[i] * (difftmp*difftmp.transpose());
      }
    cmuplus *= 1.0/(solutions._sigma*solutions._sigma);
        
    // Cmu-, Eq. (7)
    dMat cmuminus = dMat::Zero(parameters._dim,parameters._dim);
    for (int i=0;i<parameters._mu;i++)
      {
	dVec ytmp = solutions._candidates.at(parameters._lambda-i-1)._x-solutions._xmean;
	//dVec yl = (solutions._csqinv * (solutions._candidates.at(parameters._lambda-parameters._mu+i)._x-solutions._xmean)).norm() / (solutions._csqinv * ytmp).norm() * ytmp * 1.0/solutions._sigma;
	dVec yl = ytmp * 1.0/solutions._sigma; // NH says this is a good enough value.
	cmuminus += parameters._weights[i] * yl*yl.transpose();
      }
    
    // covariance update, Eq. (8)
    dMat cminusdenom = solutions._csqinv*cmuminus*solutions._csqinv;
    SelfAdjointEigenSolver<dMat> tmpesolve(cminusdenom); // XXX: computing eigenvalues, could be avoid with an upper bound.
    double cminustmp = tmpesolve.eigenvalues().maxCoeff();
    double cminusmin = parameters._alphaminusmin * (1.0-parameters._cmu)*(1.0-parameters._lambdamintarget) / cminustmp;
    double cminus = std::min(cminusmin,(1-parameters._cmu)*parameters._alphacov/8.0*(parameters._muw/(pow(parameters._dim+2.0,1.5)+2.0*parameters._muw)));
    //std::cerr << "cminus=" << cminus << " / cminusmin=" << cminusmin << " / cminustmp=" << cminustmp << std::endl;
    solutions._cov = (1-parameters._c1-parameters._cmu + cminus*parameters._alphaminusold)*solutions._cov + parameters._c1*spc + (parameters._cmu + cminus * (1.0-parameters._alphaminusold))*cmuplus - cminus * cmuminus;
    
    // sigma update, Eq. (9)
    solutions._sigma *= std::min(parameters._deltamaxsigma,exp((parameters._csigma / parameters._dsigma) * (norm_ps / parameters._chi - 1.0)));
        
    // set mean.
    solutions._xmean = xmean;
  }
  
}
