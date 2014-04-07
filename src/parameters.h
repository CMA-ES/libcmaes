
#ifndef PARAMETERS_H
#define PARAMETERS_H

namespace libcmaes
{

  class Parameters
  {
  public:
  Parameters():_dim(0),_lambda(0),_max_iter(0)
      {};
  Parameters(const int &dim, const int &lambda, const int &max_iter)
    :_dim(dim),_lambda(lambda),_max_iter(max_iter),_quiet(false)
    {};
    ~Parameters() {}

    int _dim; /**< function space dimensions. */
    int _lambda; /**< number of offsprings. */
    int _max_iter; /**< max iterations. */

    bool _quiet; /**< quiet all outputs. */
  };
  
}

#endif
