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

#ifndef SURROGATESTRATEGY_H
#define SURROGATESTRATEGY_H

#include "eo_matrix.h"
#include "cmastrategy.h"
#include <typeinfo>
#include <random>

namespace libcmaes
{
  /**
   * \brief function to train a surrogate model
   * @param candidates set of points along with objective function value
   * @param cov a possibly empty covariance matrix in order to re-scale points before training
   * @return training status
   */
  typedef std::function<int (const std::vector<Candidate>&, const dMat&)> CSurrFunc;

  /**
   * \brief function to predict from a surrogate model
   * @param candidates set of points for which value is to be predicted
   * @param cov a possibly empty covariance matrix in order to re-scale points before predicting
   * @return prediction status
   */
  typedef std::function<int (std::vector<Candidate>&, const dMat&)> SurrFunc; //TODO: a signature closer to the objective function signature ?


  /**
   * \brief Surrogate base class, to be derived in order to create strategy
   *        to be used along with CMA-ES.
   */
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class SurrogateStrategy : public TStrategy<TCovarianceUpdate,TGenoPheno>
  {
    public:
    /**
     * \brief constructor
     * @param func objective function to minimize
     * @param parameters optimization parameters
     */
    SurrogateStrategy(FitFunc &func,
		      CMAParameters<TGenoPheno> &parameters);

    protected:
    ~SurrogateStrategy();

    public:
    /**
     * \brief train a surrogate model
     * @param candidates set of points along with objective function value
     * @param cov a possibly empty covariance matrix in order to re-scale points before training
     * @return training status
     */
    int train(const std::vector<Candidate> &candidates,
	      const dMat &cov) { return _train(candidates,cov); }

    /**
     * \brief predict from a surrogate model
     * @param candidates set of points for which value is to be predicted
     * @param cov a possibly empty covariance matrix in order to re-scale points before predicting
     * @return prediction status
     */
    int predict(std::vector<Candidate> &candidates,
		const dMat &cov) { return _predict(candidates,cov); }

    /**
     * \brief compute surrogate model error (copies and sorts the test_set)
     * @param test_set the candidate points along with their objective function values for model evaluation
     * @param cov possibly empty covariance matrix in order to re-scale the points before error estimation
     * @return surrogate model error estimate
     */
    double compute_error(const std::vector<Candidate> &test_set,
			 const dMat &cov=dMat(0,0));
    
    /**
     * \brief conditionals on training, to be specialized in inherited surrogate strategies
     * @return whether to train surrogate
     */
    bool do_train() const { return true; };

    /**
     * \brief sets the training function
     * @param training function
     */
    void set_ftrain(const CSurrFunc &train) { _train = train; }

    /**
     * \brief sets the prediction function
     * @param prediction function
     */
    void set_fpredict(const SurrFunc &predict) { _predict = predict; }

    /**
     * \brief sets the size of the training set (number of points)
     * @param l size of the training set
     */
    void set_l(const int &l) { _l = l; }

    /**
     * \brief gets the size of the training set (number of points)
     * @return size of the training set
     */
    int get_l() const { return _l; }

    /**
     * \brief sets whether to exploit the surrogate model
     * @param exploit whether to exploit the surrogate model
     */
    void set_exploit(const bool &exploit) { _exploit = exploit; }

    /**
     * \brief gets the state of surrogate model exploitation
     * @return whether the surrogate model is being exploited
     */
    bool get_exploit() const { return _exploit; }

    /**
     * \brief returns the surrogate model training error
     * @return surrogate model training error
     */
    double get_train_error() const { return _train_err; }

    /**
     * \brief returns the surrogate model test error
     * @return surrogate model test error
     */
    double get_test_error() const { return _test_err; }

    /**
     * \brief sets training error
     * @param err training error
     */
    void set_train_error(const double &err) { _train_err = err; }

    /**
     * \brief sets the test error and updates the smoothed test err.
     * @param err test error
     */
    void set_test_error(const double &err);
    
    /**
     * \brief adds a point to the training set (candidate = points + objective function value)
     * @param c point to add to the training set
     */
    void add_to_training_set(const Candidate &c);

    /**
     * \brief sets the lifelength of the surrogate, i.e. the number of steps in between to training steps
     * @param nsteps surrogate lifelength, -1 for automatic determination
     */
    inline void set_nsteps(const int &nsteps)
    {
      _nsteps = nsteps;
      if (_nsteps < 0)
	_auto_nsteps = true;
    }

    /**
     * \brief resets training set and related information, useful when using algorithms with restarts
     */
    inline void reset_training_set()
    {
      _tset.clear();
      _train_err = _test_err = 0.0;
      _smooth_test_err = 0.5;
    }
    
    /**
     * \brief returns the current surrogate lifelength
     * @return current surrogate lifelength
     */
    int get_nsteps() const { return _nsteps; }
    
  protected:
    bool _exploit = true; /**< whether to exploit or test the surrogate. */
    int _l = 200; /**< number of training samples. set to floor(30*sqrt(n)) in constructor. */
    std::vector<Candidate> _tset; /**< current training set. */
    CSurrFunc _train; /**< custom training function. */
    SurrFunc _predict; /**< custom prediction function. */
    double _train_err = 0.0; /**< current surrogate training error. */
    double _test_err = 0.0; /**< current surrogate model error estimate. */
    double _smooth_test_err = 0.5; /**< smoothed test error as (1-\beta_err)*_test_err + \beta_err * new_test_err */
    double _beta_err = 0.2; /**< smoothing constant. */
    int _nsteps = 1; /**< steps in between two training phases. */
    int _auto_nsteps = false; /**< whether to automatically set the surrogate lifelength. */
    };

  /**
   * \brief Simple surrogate strategy: trains every n steps, and exploits in between,
   *        mostly as an example and for testing / debugging surrogates.
   *        This strategy overrides the ask/eval/tell functions of the base optimization strategy
   */
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
    class SimpleSurrogateStrategy : public SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>
    {
    public:
    /**
     * \brief constructor
     * @param func objective function to minimize
     * @param parameters optimization parameters
     */
    SimpleSurrogateStrategy(FitFunc &func,
			    CMAParameters<TGenoPheno> &parameters);

    ~SimpleSurrogateStrategy();

    /**
     * \brief Evaluates a set of candiates against the objective function 
     *        or the surrogate model, as needed
     *
     * Note: this function overrides the default CMAStrategy::eval
     *
     * @param candidates A matrix whose rows contain the candidates.
     * @param phenocandidates The candidates transformed into phenotype, 
     *        leave empty if no pheno transform.
     */
    void eval(const dMat &candidates,
	      const dMat &phenocandidates=dMat(0,0));

    /**
     * \brief Updates the state of the stochastic search, and prepares
     *        for the next iteration by training the surrogate model, as needed.
     *
     * Note: this function overrides the default CMAStrategy::tell
     */
    void tell();

    /**
     * \brief Finds the minimum of the objective function. It makes
     *        alternate calls to ask(), tell() and stop() until 
     *        one of the termination criteria triggers.
     * @return success or error code, as defined in opti_err.h
     * Note: the termination criteria code is held by _solutions._run_status
     */
    int optimize();
    
    /**
     * \brief estimates surrogate lifelength
     * @return estimated surrogate lifelength
     */
    int compute_lifelength();
    
    /**
     * \brief whether to train the model
     * @return whether to train the model
     */
    inline bool do_train() const
    {
      if (!SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::_exploit && (int)this->_tset.size() >= this->_l)
	return true;
      return ((this->_niter == 0 || this->_niter % this->_nsteps == 0) && (int)this->_tset.size() >= this->_l);
    }

    public:
    double _terr = 0.45; /**< error threshold for estimating optimal nsteps */
    int _nmax = 20;
    };

  /**
   * \brief ACM Surrogate strategy for CMA-ES, follows:
   *        'Surrogate-Assisted Evolutionary Algorithms', Ilya Loshchilov, PhD Thesis, Universite Paris-Sud 11, 2013.
   *        http://www.loshchilov.com/phd.html
   *        see Chapter 4.
   *
   *        Implements a single-objective strategy for CMA-ES and related algorithms, that
   *        samples ans pre-screens a larger than usual number of offsprings at each generation,
   *        rank them with a rank-based surrogate model, and consumes a small portion of offsprings
   *        with the original (supposedly expensive) objective function.
   *
   *        This strategy overrides the ask/eval/tell functions of the base optimization strategy
   */
  template<template <class U,class V> class TStrategy, class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
    class ACMSurrogateStrategy : public SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>
    {
    public:
    /**
     * \brief constructor
     * @param func objective function to minimize
     * @param parameters optimization parameters
     */
    ACMSurrogateStrategy(FitFunc &func,
			 CMAParameters<TGenoPheno> &parameters);
    
    ~ACMSurrogateStrategy();
    
    /**
     * \brief Generates a set of candidate points. Uses the pre-sampling of a larger
     *        than usual number of offprings, controled by 'lambdaprime', as needed
     * 
     * Note: this function overrides the default ESOStrategy::ask
     * 
     * @return A matrix whose rows contain the candidate points.
     */
    dMat ask();

    /**
     * \brief Evaluates a set of candiates against the objective function 
     *        or the surrogate model, as needed
     *
     * Note: this function overrides the default CMAStrategy::eval
     *
     * @param candidates A matrix whose rows contain the candidates.
     * @param phenocandidates The candidates transformed into phenotype, 
     *        leave empty if no pheno transform.
     */
    void eval(const dMat &candidates,
	      const dMat &phenocandidates=dMat(0,0));

    /**
     * \brief Updates the state of the stochastic search, and prepares
     *        for the next iteration by training the surrogate model, as needed.
     *
     * Note: this function overrides the default CMAStrategy::tell
     */
    void tell();

    int optimize();
    
    protected:
    /**
     * \brief pre-selection + candidate evaluation scheme.
     *        Called by eval, evaluates lambdaprime candidates with surrogate model,
     *        then subsample the population in order to evaluate them with the original
     *        objective function, test the surrogate model and grow the training set with 
     *        new points
     * @param candidates A matrix whose rows contain the candidates.
     */
    void pre_selection_eval(const dMat &candidates);
    
    public:
    /**
     * \brief whether to train the model
     * @return whether to train the model
     */
    inline bool do_train() const
    {
      if (!SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::_exploit && (int)this->_tset.size() >= this->_l)
	return true;
      else if (SurrogateStrategy<TStrategy,TCovarianceUpdate,TGenoPheno>::_exploit && (int)this->_tset.size() >= this->_l)
	return true;
      return ((this->_niter == 0 || this->_niter % this->_nsteps == 0) && (int)this->_tset.size() >= this->_l);
    }

    private:
    void init_rd(); // initialize random device.

    public:
    /**
     * \brief sets the number of true objective function calls per iteration
     * @param lp true objective function calls per iteration
     */
    void set_lambdaprime(const int &lp) { _lambdaprime = lp; }

    /**
     * \brief returns the number of calls to the true objective function per iteration
     * @return calls per iteration to the true objective function
     */
    int get_lambdaprime() const { return _lambdaprime; }

    /**
     * \brief sets the number of pre-screened offsprings (sampled)
     * @param number of offsprings
     */
    void set_prelambda(const int &pl) { _prelambda = pl; }

    /**
     * \brief returns the current number of pre-screened offpsrings (sampled)
     * @return number of offsprings
     */
    int get_prelambda() const { return _prelambda; }
    
    /**
     * \brief sets the standard deviation of selection sampling step 0
     * @param s standard deviation
     */
    void set_theta_sel0(const double &s) { _theta_sel0 = s; }
    
    /**
     * \brief returns the standard deviation of selection sampling step 0
     * @return standard deviation
     */
    double get_theta_sel0() const { return _theta_sel0; }

    /**
     * \brief sets the standard deviation of selection sampling step 1
     * @param s standard deviation
     */
    void set_theta_sel1(const double &s) { _theta_sel1 = s; }

    /**
     * \brief returns the standard deviation of selection sampling step 0
     * @return standard deviation
     */
    double get_theta_sel1() const { return _theta_sel1; }

    protected:
    double _prelambda = 500; /**< number of pre-screened offsprings. */
    double _theta_sel0 = 0.4;  /**< standard deviation of selection sampling step 0. */
    double _theta_sel1 = 0.8;  /**< standard deviation of selection sampling step 1. */
    int _lambdaprime; /**< true objective function calls per iteration. */

    private:
    // random numbers for selection sampling
    std::random_device _rd;
    std::normal_distribution<double> _norm_sel0;
    std::normal_distribution<double> _norm_sel1;
    std::mt19937 _gen0;
    std::mt19937 _gen1;
  };
}

#endif
