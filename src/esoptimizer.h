
#ifndef ESOPTIMIZER_H
#define ESOPTIMIZER_H

#include <functional>
#include <chrono>
#include "parameters.h"
#include "cmastrategy.h"

namespace libcmaes
{
  
  template <class TESOStrategy,class TParameters>
    class ESOptimizer : public TESOStrategy
    {
    public:
      ESOptimizer(FitFunc &func,
		  TParameters &parameters)
	:TESOStrategy(func,parameters),_elapsed_optimize_ms(0)
	{
	}
      
      ~ESOptimizer() {};

      bool optimize()
      {
	std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
	bool opt = TESOStrategy::optimize();
	std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	_elapsed_optimize_ms = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	return opt;
      }

      int _elapsed_optimize_ms;
    };

}

#endif
