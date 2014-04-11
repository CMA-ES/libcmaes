
#ifndef IPOPCMASTRATEGY_H
#define IPOPCMASTRATEGY_H

#include "cmastrategy.h"
#include "cmaparameters.h"

namespace libcmaes
{
  class IPOPCMAStrategy : public CMAStrategy<CovarianceUpdate>
  {
  public:
    IPOPCMAStrategy(FitFunc &func,
		    CMAParameters &parameters);
    ~IPOPCMAStrategy();

    void tell();

    int optimize();
    
  };
}

#endif
