#include "esoptimizer.h"
#include "cmastrategy.h"
#include <glog/logging.h>

using namespace libcmaes;

FitFunc cigtab = [](const double *x, const int N)
{
  int i;
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)
    sum += x[i]*x[i];
  return sum;
};

int main(int argc, char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr=1;
  google::SetLogDestination(google::INFO,"");
  //FLAGS_log_prefix=false;

  int dim = 5;
  int lambda = 10;
  CMAParameters cmaparams(dim,lambda);
  ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters> cmaes(cigtab,cmaparams);
  cmaes.optimize();
}
