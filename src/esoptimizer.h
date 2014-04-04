
#ifndef ESOPTIMIZER_H
#define ESOPTIMIZER_H

#include <functional>
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
	:TESOStrategy(func,parameters)
	{
	}
      /*ESOptimizer(FitFunc &func,
		  const int &dim,
		const int &lambda) //deprecated
	:TESOStrategy(func,dim,lambda)
      {
      };*/
      
      ~ESOptimizer() {};

      bool optimize() { return TESOStrategy::optimize(); };
    };

}

#endif
