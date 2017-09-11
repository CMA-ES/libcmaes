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

#ifndef BIPOPCMASTRATEGY_H
#define BIPOPCMASTRATEGY_H

#include "ipopcmastrategy.h"
#include <random>

namespace libcmaes
{
  /**
   * \brief Implementation of the BIPOP flavor of CMA-ES, with restarts that
   *        control the population of offsprings used in the update of the 
   *        distribution parameters in order to alternate between local and 
   *        global searches for the objective.
   */
  template <class TCovarianceUpdate, class TGenoPheno>
    class CMAES_EXPORT BIPOPCMAStrategy : public IPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>
  {
  public:
    /**
     * \brief constructor.
     * @param func objective function to minimize
     * @param parameters stochastic search parameters
     */
    BIPOPCMAStrategy(FitFunc &func,
		     CMAParameters<TGenoPheno> &parameters);

    /**
     * \brief constructor.
     * @param func objective function to minimize
     * @param parameters stochastic search parameters
     * @param solutions solution to start search from
     */
    BIPOPCMAStrategy(FitFunc &func,
		     CMAParameters<TGenoPheno> &parameters,
		     const CMASolutions &solutions);
    
    ~BIPOPCMAStrategy();

    /**
     * \brief Updates the covariance matrix and prepares for the next iteration.
     */
    void tell();

    /**
     * \brief Finds the minimum of the objective function. It makes
     *        alternate calls to ask(), tell() and stop() until 
     *        one of the termination criteria triggers.
     * @return success or error code, as defined in opti_err.h
     * Note: the termination criteria code is held by _solutions._run_status
     */
    int optimize(const EvalFunc &evalf,const AskFunc &askf,const TellFunc &tellf);

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
    int optimize()
    {
      return optimize(std::bind(&BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::eval,this,std::placeholders::_1,std::placeholders::_2),
		      std::bind(&BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::ask,this),
		      std::bind(&BIPOPCMAStrategy<TCovarianceUpdate,TGenoPheno>::tell,this));
    }
    
  protected:
    void r1();
    void r2();

  private:
    std::mt19937 _gen;
    std::uniform_real_distribution<> _unif;
    double _lambda_def;
    double _lambda_l;
    double _sigma_init; // to save the original value
    double _max_fevals; // to save the original value
  };
}

#endif
