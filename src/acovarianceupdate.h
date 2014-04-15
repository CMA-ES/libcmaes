
#ifndef ACOVARIANCEUPDATE_H
#define ACOVARIANCEUPDATE_H

#include "cmaparameters.h"
#include "cmasolutions.h"
#include "eigenmvn.h"

namespace libcmaes
{

  //TODO: paper reference.
  
  class ACovarianceUpdate
  {
  public:
    static void update(const CMAParameters &parameters,
		       EigenMultivariateNormal<double> &esolver,
		       CMASolutions &solutions);
  };
  
}

#endif
