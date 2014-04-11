
#ifndef ESOSTRATEGY_H
#define ESOSTRATEGY_H

#include "eo_matrix.h" // to include Eigen everywhere.
#include "candidate.h"

namespace libcmaes
{
  typedef std::function<double (const double*, const int &n)> FitFunc;

  template<class TParameters,class TSolutions,class TStopCriteria>
    class ESOStrategy
  {
  public:
    ESOStrategy(FitFunc &func,
		TParameters &parameters);

  protected:
    ~ESOStrategy();

  public:
    /**
     * \brief generates nsols new candidate solutions, sampled from a 
     *        multivariate normal distribution.
     * TODO: move to CMAStrategy as specific to CMA-ES?
     */
    dMat ask();

    void eval(const dMat &candidates);

    void tell();
    
    bool stop();

    int optimize();

    // deprecated.
    Candidate best_solution() const;
    
    FitFunc _func;
    int _nevals;  /**< number of function evaluations. */
    int _niter;  /**< number of iterations. */
    TSolutions _solutions;
    TParameters _parameters;
  };
  
}

#endif
