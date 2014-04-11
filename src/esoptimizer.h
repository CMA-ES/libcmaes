
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
  /**
   * \brief an optimizer main class. 
   */
  template <class TESOStrategy,class TParameters>
    class ESOptimizer : public TESOStrategy
    {
    public:
      /**
       * \brief constructor
       * @param func function to minimize
       * @param parameters optimization parameters
       */
      ESOptimizer(FitFunc &func,
		  TParameters &parameters)
	:TESOStrategy(func,parameters)
	{
	}
      
      ~ESOptimizer() {};

      /**
       * \brief finds the minimum of a function, by calling on the underlying
       *        procedure of the EOSOptimizer object, like a variety of flavor of CMA-ES.
       */
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
