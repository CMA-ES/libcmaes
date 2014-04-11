
#ifndef ESOPTIMIZER_H
#define ESOPTIMIZER_H

#include <functional>
#include <chrono>
#include "parameters.h"
#include "cmastrategy.h"

/* algorithms */
enum {
  /* vanilla version of CMA-ES. */
  CMAES_DEFAULT = 0,
  /* IPOP-CMA-ES. */
  IPOP_CMAES = 1
};

namespace libcmaes
{
  
  template <class TESOStrategy,class TParameters>
    class ESOptimizer : public TESOStrategy
    {
    public:
      ESOptimizer(FitFunc &func,
		  TParameters &parameters)
	:TESOStrategy(func,parameters)
	{
	}
      
      ~ESOptimizer() {};

      
      
      int optimize()
      {
	std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
	bool opt = TESOStrategy::optimize();
	std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	TESOStrategy::_solutions._elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	return opt;
      }
    };

}

#endif
