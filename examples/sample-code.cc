#include "cmaes.h"
#include <iostream>

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

int main(int argc, char *argv[])
{
  int dim = 10; // problem dimensions.
  //int lambda = 100; // offsprings at each generation.
  //CMAParameters cmaparams(dim,lambda);
  CMAParameters<> cmaparams(dim);
  //cmaparams._algo = BIPOP_CMAES;
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  std::cout << "best solution: " << cmasols << std::endl;
  std::cout << "optimization took " << cmasols._elapsed_time / 1000.0 << " seconds\n";
  return cmasols._run_status;
}
