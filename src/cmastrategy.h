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

#ifndef CMASTRATEGY_H
#define CMASTRATEGY_H

#include "esostrategy.h"
#include "cmaparameters.h"
#include "cmasolutions.h"
#include "cmastopcriteria.h"
#include "covarianceupdate.h"
#include "acovarianceupdate.h"
#include "vdcmaupdate.h"
#include "eigenmvn.h"
#include <fstream>

namespace libcmaes
{

  /**
   * \brief This is an implementation of CMA-ES. It uses the reference algorithm
   *        and termination criteria of the following paper:
   * Hansen, N. (2009). Benchmarking a BI-Population CMA-ES on the BBOB-2009 Function Testbed. Workshop Proceedings of the GECCO Genetic and Evolutionary Computation Conference, ACM, pp. 2389-2395
   * See https://www.lri.fr/~hansen/publications.html for more information.
   */
  template <class TCovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
    class CMAStrategy : public ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >
    {
    public:
      /**
       * \brief dummy constructor
       */
      CMAStrategy();
    
      /**
       * \brief constructor.
       * @param func objective function to minimize
       * @param parameters stochastic search parameters
       */
      CMAStrategy(FitFunc &func,
		  CMAParameters<TGenoPheno> &parameters);

      /**
       * \brief constructor for starting from an existing solution.
       * @param func objective function to minimize
       * @param parameters stochastic search parameters
       * @param cmasols solution object to start from
       */
      CMAStrategy(FitFunc &func,
		  CMAParameters<TGenoPheno> &parameters,
		  const CMASolutions &cmasols);
    
      ~CMAStrategy();

      /**
       * \brief generates nsols new candidate solutions, sampled from a 
       *        multivariate normal distribution.
       * return A matrix whose rows contain the candidate points.
       */
      dMat ask();

      /**
       * \brief Updates the covariance matrix and prepares for the next iteration.
       */
      void tell();

      /**
       * \brief Stops search on a set of termination criterias, see reference paper.
       * @return true if search must stop, false otherwise.
       */
      bool stop();

      /**
       * \brief Finds the minimum of the objective function. It makes
       *        alternate calls to ask(), tell() and stop() until 
       *        one of the termination criteria triggers.
       * @param evalf custom eval function
       * @param askf custom ask function
       * @param tellf custom tell function
       * @return success or error code, as defined in opti_err.h
       * Note: the termination criteria code is held by _solutions._run_status
       */
    int optimize(const EvalFunc &evalf, const AskFunc &askf, const TellFunc &tellf);

      /**
       * \brief Finds the minimum of the objective function. It makes
       *        alternate calls to ask(), tell() and stop() until 
       *        one of the termination criteria triggers.
       * @return success or error code, as defined in opti_err.h
       * Note: the termination criteria code is held by _solutions._run_status
       */
    int optimize()
    {
      return optimize(std::bind(&ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno>>::eval,this,std::placeholders::_1,std::placeholders::_2),
		      std::bind(&CMAStrategy<TCovarianceUpdate,TGenoPheno>::ask,this),
		      std::bind(&CMAStrategy<TCovarianceUpdate,TGenoPheno>::tell,this));
    }

      /**
       * \brief Stream the internal state of the search into an output file, 
       *        as defined in the _parameters object.
       */
      void plot();
    
    protected:
      Eigen::EigenMultivariateNormal<double> _esolver;  /**< multivariate normal distribution sampler, and eigendecomposition solver. */
      CMAStopCriteria<TGenoPheno> _stopcriteria; /**< holds the set of termination criteria, see reference paper. */
      std::ofstream *_fplotstream = nullptr; /**< plotting file stream, not in parameters because of copy-constructor hell. */
    
    public:
    static ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> _defaultPFunc; /**< the default progress function. */
    static PlotFunc<CMAParameters<TGenoPheno>,CMASolutions> _defaultFPFunc; /**< the default plot to file function. */
    };
  
}

#endif
