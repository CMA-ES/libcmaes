
#ifndef CMASTRATEGY_H
#define CMASTRATEGY_H

#include "esostrategy.h"
#include "covarianceupdate.h"
#include "cmaparameters.h"
#include "cmasolutions.h"
#include "cmastopcriteria.h"
#include "eigenmvn.h"

namespace libcmaes
{
  template <class TCovarianceUpdate>
    class CMAStrategy : public ESOStrategy<CMAParameters,CMASolutions,CMAStopCriteria>
    {
    public:
      CMAStrategy(FitFunc &func,
		  CMAParameters &parameters);
      ~CMAStrategy();

      dMat ask();

      void tell();

      bool stop();

      bool optimize();
      
    private:
      EigenMultivariateNormal<double> _esolver;
      CMAStopCriteria _stopcriteria;
    };
  
}

#endif
