#include "esoptimizer.h"
#include "cmastrategy.h"
#include "llogging.h"

using namespace libcmaes;

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
};

GradFunc grad_fsphere = [](const double *x, const int N)
{
  dVec grad(N);
  for (int i=0;i<N;i++)
    grad(i) = 2.0*x[i];
  return grad;
};

FitFunc elli = [](const double *x, const int N)
{
  if (N == 1)
    return x[0] * x[0];
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += exp(log(1e3)*2.0*static_cast<double>(i)/static_cast<double>((N-1))) * x[i]*x[i];
  return val;
};

FitFunc rosenbrock = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N-1;i++)
    {
      val += 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
    }
  return val;
};

int main(int argc, char *argv[])
{
  int samples = 100;
  std::string sep = "\t";
  std::cout << "lambda\tEDM\tEDM/f(m)\n";
  std::vector<int> lambdas = {5,10,20,40,80,160,320,640,1280};
  int dim = 10;
  std::vector<double> x0(dim,1.0);
  double sigma = 100.0;
  int s = 0;
  while (s < samples)
    {
      for (size_t i=0;i<lambdas.size();i++)
	{
	  double lambda = lambdas.at(i);
	  CMAParameters<> cmaparams(x0,sigma,lambda);
	  cmaparams.set_quiet(true);
	  //cmaparams.set_ftarget(1e-8);
	  //cmaparams.set_algo(aCMAES);
	  ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<>> cmaes(fsphere,cmaparams);
	  //cmaes.set_gradient_func(grad_fsphere);
	  cmaes.optimize();
	  std::cout << lambda << sep << cmaes.get_solutions().edm() << sep << cmaes.get_solutions().edm() /  cmaes.fitfunc(cmaes.get_solutions().xmean().data(),dim) << std::endl;
	}
      ++s;
    }
}
