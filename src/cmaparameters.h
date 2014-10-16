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

#ifndef CMAPARAMETERS_H
#define CMAPARAMETERS_H

#include "parameters.h"
#include "eo_matrix.h"
#include <float.h>

namespace libcmaes
{
  /**
   * \brief Parameters for various flavors of the CMA-ES algorithm.
   */
  template <class TGenoPheno=GenoPheno<NoBoundStrategy> >
  class CMAParameters : public Parameters<TGenoPheno>
    {
      friend class CMASolutions;
      template <class U, class V> friend class CMAStrategy;
      template <class U, class V, class W> friend class ESOStrategy;
      template <class U> friend class CMAStopCriteria;
      template <class U, class V> friend class IPOPCMAStrategy;
      template <class U, class V> friend class BIPOPCMAStrategy;
      friend class CovarianceUpdate;
      friend class ACovarianceUpdate;
      
    public:
      CMAParameters() {} //TODO: var init even if this constructor is not supposed to be used for now.
      
      /**
       * \brief Constructor.
       * @param dim problem dimensions
       * @param x0 initial search point
       * @param sigma initial distribution step size (positive, otherwise automatically set)
       * @param lambda number of offsprings sampled at each step
       * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
       * @param gp genotype / phenotype object
       * @param sep whether to use sep-CMA-ES, using diagonal covariance matrix (modifies covariance default learning rate)
       */
      CMAParameters(const int &dim,
		    const double *x0,
		    const double &sigma,
		    const int &lambda=-1,
		    const uint64_t &seed=0,
		    const TGenoPheno &gp=GenoPheno<NoBoundStrategy>());
      
      /**
       * \brief Constructor.
       * @param x0 initial search point as vector of problem dimension
       * @param sigma initial distribution step size (positive, otherwise automatically set)
       * @param lambda number of offsprings sampled at each step
       * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
       * @param gp genotype / phenotype object
       * @param sep whether to use sep-CMA-ES, using diagonal covariance matrix (modifies covariance default learning rate)
       */
      CMAParameters(const std::vector<double> &x0,
		    const double &sigma,
		    const int &lambda=-1,
		    const uint64_t &seed=0,
		    const TGenoPheno &gp=GenoPheno<NoBoundStrategy>());
      
      /**
       * \brief Constructor.
       * @param x0 initial search point as vector of problem dimension
       * @param sigma vector of initial distribution step sizes (positive, otherwise automatically set)
       * @param lambda number of offsprings sampled at each step
       * @param seed initial random seed, useful for reproducing results (if unspecified, automatically generated from current time)
       * @param gp genotype / phenotype object
       * @param sep whether to use sep-CMA-ES, using diagonal covariance matrix (modifies covariance default learning rate)
       */
      CMAParameters(const std::vector<double> &x0,
		    const std::vector<double> &sigma,
		    const int &lambda=-1,
		    const std::vector<double> &lbounds=std::vector<double>(),
		    const std::vector<double> &ubounds=std::vector<double>(),
		    const uint64_t &seed=0);
      
      ~CMAParameters();
      
      /**
       * \brief initialize required parameters based on dim, lambda, x0 and sigma.
       */
      void initialize_parameters();
      
      /**
       * \brief adapt parameters for noisy objective function.
       */
      void set_noisy();
      
      /**
       * \brief sets the optimization algorithm.
       * @param algo as string from cmaes,ipop,bipop,acmaes,aipop,abipop,sepcmaes,sepipop,sepbipop,sepacmaes,sepaipop,sepabipop
       */
      void set_str_algo(const std::string &algo)
      {
	std::map<std::string,int>::const_iterator mit;
	if ((mit = Parameters<TGenoPheno>::_algos.find(algo))!=Parameters<TGenoPheno>::_algos.end())
	  Parameters<TGenoPheno>::_algo = (*mit).second;
	else LOG(ERROR) << "unknown algorithm " << algo << std::endl;
	if (algo.find("sep")!=std::string::npos)
	  set_sep();
      }
      
      /**
       * \brief fix parameters for sep-CMA-ES, using only the diagonal of covariance matrix.
       */
      void set_sep();

      /**
       * \brief whether algorithm leverages separability.
       * @return separability status
       */
      bool is_sep() const { return _sep; }
      
      /**
       * \brief turns stopping criteria MaxIter that automatically stops optimization after a 
       *        number of steps on or off.
       * @param b true or false for turning criteria on or off (on is default in constructor).
       */
      inline void set_automaxiter(const bool &b) { _has_max_iter = b; }
      
      /**
       * \brief freezes a parameter to a given value during optimization.
       *        Adapts some generic parameters as well.
       * @param index dimension index of the parameter to be frozen
       * @param value frozen value of the parameter
       */
      void set_fixed_p(const int &index, const double &value);
      
      /**
       * \brief unfreezes a parameter.
       * @param index dimenion index of the parameter to unfreeze
       */
      void unset_fixed_p(const int &index);
      
      /**
       * \brief sets the maximum number of restarts (applies to IPOP and BIPOP).
       * @param nrestarts maximum number of restarts
       */
      inline void set_restarts(const int &nrestarts) { _nrestarts = nrestarts; }

      /**
       * \brief get the number of restarts (applies to IPOP and BIPOP).
       * @return number of restarts
       */
      inline int get_restarts() const { return _nrestarts; }

      /**
       * \brief sets the lazy update (i.e. updates the eigenvalues every few steps).
       * @param lz whether to activate the lazy update
       */
      inline void set_lazy_update(const bool &lz) { _lazy_update = lz; }

      /**
       * \brief get lazy update status.
       * @param whether lazy update is activated
       */
      inline bool get_lazy_update() { return _lazy_update; }

      /**
       * \brief sets initial elitist scheme: restart if best encountered solution is not
       *        the final solution and reinjects the best solution until the population
       *        has better fitness, in its majority
       * @param e whether to activate the initial elitist scheme
       */
      inline void set_elitist(const bool &e) { _elitist = e; }

      /**
       * \brief sets the optima hoppping scheme: restart while best encountered solution is not
       *        the final solution and reinects the best solution during each run until the population
       *        has better fitness, in its majority
       * @param oh whether to activate the optimal hopping scheme
       */
      inline void set_optima_hopping(const bool &oh) { _optima_hopping = oh; }
      
    private:
      int _mu; /**< number of candidate solutions used to update the distribution parameters. */
      dVec _weights; /**< offsprings weighting scheme. */
      double _csigma; /**< cumulation constant for step size. */
      double _c1; /**< covariance matrix learning rate for the rank one update using pc. */
      double _cmu; /**< covariance matrix learning reate for the rank mu update. */
      double _cc; /**< cumulation constant for pc. */
      double _muw; /**< \sum^\mu _weights .*/
      double _dsigma; /**< step size damping factor. */
      
      // computed once at init for speeding up operations.
      double _fact_ps;
      double _fact_pc;
      double _chi; /**< norm of N(0,I) */
      
      double _sigma_init; /**< initial sigma value. */
      
      int _nrestarts = 9; /**< maximum number of restart, when applicable. */
      bool _lazy_update; /**< covariance lazy update. */
      double _lazy_value; /**< reference trigger for lazy update. */
      
      // active cma.
      double _cm; /**< learning rate for the mean. */
      double _alphacov; /**< = 2 (active CMA only) */
      double _alphaminusold; /**< in [0,1] (active CMA only) */
      double _deltamaxsigma; /**< infinite (active CMA only) */
      double _lambdamintarget; /**< = 0.66 (active CMA only) */
      double _alphaminusmin; /**< = 1 (active CMA only) */
      
      // sep cma (diagonal cov).
      bool _sep = false; /**< whether to use diagonal covariance matrix. */
      
      // stopping criteria.
      bool _has_max_iter = true; /**< MaxIter criteria: automatically stop running after 100+50*((D+2)^2)/lambda iterations. */

      bool _elitist = false; /**< activate the restart from and re-injection of the best seen solution if not the final one. */
      bool _optima_hopping = false; /**< activates a form of elistism that always restart from the best seen solution. */
    };

  template<class TGenoPheno>
    std::map<std::string,int> Parameters<TGenoPheno>::_algos = {{"cmaes",0},{"ipop",1},{"bipop",2},{"acmaes",3},{"aipop",4},{"abipop",5},{"sepcmaes",6},{"sepipop",7},{"sepbipop",8},{"sepacmaes",9},{"sepipop",10},{"sepbipop",11}};
}

#endif
