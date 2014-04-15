
#ifndef COVARIANCEUPDATE_H
#define COVARIANCEUPDATE_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include "eigenmvn.h"

namespace libcmaes
{

  class CovarianceUpdate
  {
  public:
    static void update(const CMAParameters &parameters,
		       EigenMultivariateNormal<double> &esolver,
		       CMASolutions &solutions);
  };
  
}

#endif
