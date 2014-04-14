
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <time.h>
#include <limits>

namespace libcmaes
{

  class Parameters
  {
  public:
  Parameters():_dim(0),_lambda(0),_max_iter(0)
      {}
  Parameters(const int &dim, const int &lambda, const int &max_iter,
	     const int &max_fevals=-1,
	     const double &x0=std::numeric_limits<double>::min(), const std::string &fplot="",
	     const uint64_t &seed=0)
    :_dim(dim),_lambda(lambda),_max_iter(max_iter),_max_fevals(max_fevals),
      _quiet(false),_fplot(fplot),_x0(x0),_seed(seed),_algo(0)
    {
      if (_seed == 0) // seed is not forced.
	_seed = static_cast<uint64_t>(time(NULL));
    }
  ~Parameters()
    {
    }

    int _dim; /**< function space dimensions. */
    int _lambda; /**< number of offsprings. */
    int _max_iter; /**< max iterations. */
    int _max_fevals; /**< max budget as number of function evaluations. */
    
    bool _quiet; /**< quiet all outputs. */
    std::string _fplot; /**< plotting file, if specified. */
    double _x0; /**< initial mean vector value for all components. */
    
    uint64_t _seed; /**< seed for random generator. */
    int _algo; /**< selected algorithm. */
  };
  
}

#endif
