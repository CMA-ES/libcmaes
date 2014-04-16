
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
  class ACovarianceUpdate
  {
  public:
    static void update(const CMAParameters &parameters,
		       EigenMultivariateNormal<double> &esolver,
		       CMASolutions &solutions);
  };
  
}

#endif
