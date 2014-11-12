/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ESOPTIMIZER_H
#define ESOPTIMIZER_H

#include <functional>
#include <chrono>
#include "parameters.h"
#include "esostrategy.h"
#include "cmasolutions.h"

/* algorithms */
enum {
  /* vanilla version of CMA-ES. */
  CMAES_DEFAULT = 0,
  /* IPOP-CMA-ES */
  IPOP_CMAES = 1,
  /* BIPOP-CMA-ES */
  BIPOP_CMAES = 2,
  /* Active CMA-ES */
  aCMAES = 3,
  /* Active IPOP-CMA-ES */
  aIPOP_CMAES = 4,
  /* Active BIPOP-CMA-ES */
  aBIPOP_CMAES = 5,
  /* sep-CMA-ES */
  sepCMAES = 6,
  /* sep-IPOP-CMA-ES */
  sepIPOP_CMAES = 7,
  /* sep-BIPOP-CMA-ES */
  sepBIPOP_CMAES = 8,
  /* Active sep-CMA-ES */
  sepaCMAES = 9,
  /* Active sep-IPOP-CMA-ES */
  sepaIPOP_CMAES = 10,
  /* Active sep-BIPOP-CMA-ES */
  sepaBIPOP_CMAES = 11,
  /* VD-CMA-ES */
  VD_CMAES = 12,
  /* VD-IPOP-CMA-ES */
  VD_IPOP_CMAES = 13,
  /* VD-BIPOP-CMA-ES */
  VD_BIPOP_CMAES = 14
};

namespace libcmaes
{
  /**
   * \brief an optimizer main class. 
   */
  template <class TESOStrategy,class TParameters,class TSolutions=CMASolutions>
    class ESOptimizer : public TESOStrategy
    {
    public:
      /**
       * \brief dummy constructor
       */
      ESOptimizer()
       :TESOStrategy()
      {
      }
    
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
      
      /**
       * \brief constructor for starting from an existing solution
       * @param func function to minimize
       * @param parameters optimization parameters
       * @param solution solution to start from
       */
      ESOptimizer(FitFunc &func,
		  TParameters &parameters,
		  const TSolutions &solution)
	:TESOStrategy(func,parameters,solution)
	{
	}
      
      ~ESOptimizer() {}

      /**
       * \brief finds the minimum of a function, by calling on the underlying
       *        procedure of the EOSOptimizer object, like a variety of flavor of CMA-ES.
       */
      int optimize()
      {
	std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
	int opt = TESOStrategy::optimize();
	std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	TESOStrategy::_solutions._elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	return opt;
      }
    };

}

#endif
