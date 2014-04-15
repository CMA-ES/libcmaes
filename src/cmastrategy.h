
#ifndef CMASTRATEGY_H
#define CMASTRATEGY_H

#include "esostrategy.h"
#include "cmaparameters.h"
//#include "acmaparameters.h"
#include "cmasolutions.h"
#include "cmastopcriteria.h"
#include "covarianceupdate.h"
#include "acovarianceupdate.h"
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
  template <class TCovarianceUpdate>
    class CMAStrategy : public ESOStrategy<CMAParameters,CMASolutions,CMAStopCriteria>
    {
    public:
      CMAStrategy(FitFunc &func,
		  CMAParameters &parameters);
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
       *        alternative calls to ask(), tell() and stop() until 
       *        one of the termination criteria triggers.
       * @return success or error code, as defined in opti_err.h
       * Note: the termination criteria code is held by _solutions._run_status
       */
      int optimize();

      /**
       * \brief Stream the internal state of the search into an output file, 
       *        as defined in the _parameters object.
       */
      void plot();

      /**
       * \brief Sets the possibly custom progress function,
       *        that is called in between every search step, and gives an outside
       *        user a simple way to witness progress of the algorithm, as well as
       *        to add custom termination criteria.
       * @param pfunc a progress function
       */
      void set_progress_func(ProgressFunc<CMAParameters,CMASolutions> &pfunc) { _defaultPFunc = pfunc; }
      
    private:
      EigenMultivariateNormal<double> _esolver;  /**< multivariate normal distribution sampler, and eigendecomposition solver. */
      CMAStopCriteria _stopcriteria; /**< holds the set of termination criteria, see reference paper. */
      std::ofstream _fplotstream; /**< plotting file stream, not in parameters because of copy-constructor hell. */

    public:
      static ProgressFunc<CMAParameters,CMASolutions> _defaultPFunc; /**< the default progress function. */
    };
  
}

#endif
