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

#include "libcmaes_config.h"
#include "esostrategy.h"
#include "cmaparameters.h" // in order to pre-instanciate template into library.
#include "cmasolutions.h"
#include "cmastopcriteria.h"
#include <iostream>
#include <glog/logging.h>

#ifdef HAVE_DEBUG
#include <chrono>
#endif

namespace libcmaes
{
  template<class TParameters,class TSolutions,class TStopCriteria>
  ESOStrategy<TParameters,TSolutions,TStopCriteria>::ESOStrategy(FitFunc &func,
						   TParameters &parameters)
    :_func(func),_nevals(0),_niter(0),_parameters(parameters)
  {
    _pfunc = [](const TParameters&,const TSolutions&){return 0;}; // high level progress function does do anything.
    _solutions = TSolutions(_parameters);
  }
  
  template<class TParameters,class TSolutions,class TStopCriteria>
  ESOStrategy<TParameters,TSolutions,TStopCriteria>::~ESOStrategy()
  {
  }
  
  template<class TParameters,class TSolutions,class TStopCriteria>
  void ESOStrategy<TParameters,TSolutions,TStopCriteria>::eval(const dMat &candidates)
  {
#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
#endif
    
    // one candidate per row.
#pragma omp parallel for if (candidates.cols() >= 100)
    for (int r=0;r<candidates.cols();r++)
      {
	_solutions._candidates.at(r)._x = candidates.col(r);
	_solutions._candidates.at(r)._fvalue = _func(_solutions._candidates.at(r)._x.data(),candidates.rows());
	
	//std::cerr << "candidate x: " << _solutions._candidates.at(r)._x.transpose() << std::endl;
      }
    _nevals += candidates.cols();
    _solutions._nevals = _nevals;

#ifdef HAVE_DEBUG
    std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
    _solutions._elapsed_eval = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
#endif
  }

  template<class TParameters,class TSolutions,class TStopCriteria>
  Candidate ESOStrategy<TParameters,TSolutions,TStopCriteria>::best_solution() const
  {
    return _solutions.best_candidate();
  }
  
  template class ESOStrategy<CMAParameters<GenoPheno<NoBoundStrategy>>,CMASolutions,CMAStopCriteria<GenoPheno<NoBoundStrategy>> >;
  template class ESOStrategy<CMAParameters<GenoPheno<pwqBoundStrategy>>,CMASolutions,CMAStopCriteria<GenoPheno<pwqBoundStrategy>> >;
}
