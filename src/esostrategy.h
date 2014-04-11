
#ifndef ESOSTRATEGY_H
#define ESOSTRATEGY_H

#include "eo_matrix.h" // to include Eigen everywhere.
#include "candidate.h"

namespace libcmaes
{
  typedef std::function<double (const double*, const int &n)> FitFunc;

  template<class TParameters,class TSolutions>
    using ProgressFunc = std::function<int (const TParameters, const TSolutions)>; // template aliasing.

  /**
   * \brief Main class describing an evolutionary optimization strategy.
   *        Every algorithm in libcmaes descends from this class, and bring
   *        its functionalities to an ESOptimizer object.
   */
  template<class TParameters,class TSolutions,class TStopCriteria>
    class ESOStrategy
  {
  public:
    /**
     * \brief constructor
     * @param func function to minimize
     * @param parameters optimization parameters
     */
    ESOStrategy(FitFunc &func,
		TParameters &parameters);

  protected:
    ~ESOStrategy();

  public:
    /**
     * \brief Generates a set of candidate points.
     * @return A matrix whose rows contain the candidate points.
     */
    dMat ask();

    /**
     * \brief Evaluates a set of candiates against the objective function.
     *        The procedure is multithreaded and stores both the candidates
     *        and their f-value into the _solutions object that bears the 
     *        current set of potential solutions to the optimization problem.
     * @param candidates A matrix whose rows contain the candidates.
     */
    void eval(const dMat &candidates);

    /**
     * \brief Updates the state of the stochastic search, and prepares
     *        for the next iteration.
     */
    void tell();

    /**
     * \brief Decides whether to stop the search for solutions.
     * @return true if search must stop, false otherwise.
     */
    bool stop();

    /**
     * \brief Finds the minimum of the objective function. It makes
     *        alternative calls to ask(), tell() and stop() until 
     *        one of the termination criteria triggers.
     * @return success or error code, as defined in opti_err.h
     */
    int optimize();

    // deprecated.
    Candidate best_solution() const;
    
    FitFunc _func; /**< the objective function. */
    int _nevals;  /**< number of function evaluations. */
    int _niter;  /**< number of iterations. */
    TSolutions _solutions; /**< holder of the current set of solutions and the dynamic elemenst of the search state in general. */
    TParameters _parameters; /**< the optimizer's set of static parameters, from inputs or internal. */
    ProgressFunc<TParameters,TSolutions> _pfunc; /**< possibly custom progress function. */
  };
  
}

#endif
