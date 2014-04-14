
#ifndef IPOPCMASTRATEGY_H
#define IPOPCMASTRATEGY_H

#include "cmastrategy.h"

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

  protected:
    void lambda_inc();
    void reset_search_state();
    void capture_best_solution(CMASolutions &best_run);
  };
}

#endif
