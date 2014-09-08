## libcmaes
libcmaes is a multithreaded C++ implementation of the CMA-ES algorithm for stochastic optimization of nonlinear 'blackbox' functions. The implemented algorithms have a wide range of applications in various disciplines, ranging from pure function minimization, optimization in industrial and scientific applications, to the solving of reinforcement and machine learning problems.

Over the past decade, both the original CMA-ES and its improved flavors have proven very effective in optimizing functions when no gradient is available. Typically, the algorithm does find the minimum value of an objective function in a minimal number of function calls, compared to other methods. For a full report of recent results, see (3).

CMA-ES is mostly the work of Nikolaus Hansen (4) and a few others. Other implementations can be found in (5).

Main functionalities:
At the moment, the library implements a vanilla version of CMA-ES (1).
Current features include:

- high-level API for simple use in external applications;
- implements several flavors of CMA-ES, IPOP-CMA-ES, BIPOP-CMA-ES, active CMA-ES, active IPOP and BIPOP restart strategies, sep-CMA-ES (linear time & space complexity) along with support for IPOP and BIPOP flavors as well;
- some operations benefit from multicores;
- support for objective function gradient, when available;
- a control exe in the command line for running the algorithm over a range of classical single-objective optimization problems.


Documentation:

- Full documentation is available from https://github.com/beniz/libcmaes/wiki
- API documentation is available from http://beniz.github.io/libcmaes/doc/html/index.html

Dependencies:

- [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for all matrix operations;
- [glog](https://code.google.com/p/google-glog/) for logging events and debug (optional);
- [gflags](https://code.google.com/p/gflags/) for command line parsing (optional);
- [gtest](https://code.google.com/p/googletest/) for unit testing (optional).

Implementation:
The library makes use of C++ policy design for modularity, performance and putting the maximum burden onto the compile-time checks. The implementation closely follows the algorithms described in (2) and (6).

### Authors
libcmaes is designed and implemented by Emmanuel Benazera on behalf of Inria Saclay / Research group TAO / LAL Appstats.

### Build
Below are instruction for Linux systems, for building on Mac, see https://github.com/beniz/libcmaes/wiki/Building-libcmaes-on-Mac-OSX

Beware of dependencies, typically on Debian/Ubuntu Linux, do:

```
sudo apt-get install autoconf automake libtool libgoogle-glog-dev libgflags-dev libeigen3-dev
```

For compiling with basic options enabled:
```
./autogen.sh
./configure
make
```

### Run examples
```
cd tests
./test_functions --dim 30 --lambda 100 --max_iter 120 --fname fsphere
```
to minimize the sphere function in 30D with 100 offsprings per generation,
```
./test_functions --dim 20 --lambda 100 --max_iter 1000 --fname rosenbrock
```
to minimize the Rosenbrock function in 20D with 100 offsprings. To see available function, do
```
./test_functions --list
```
to plot results, use a file output and then the included Gnuplot script
```
./test_functions --fname rastrigin --dim 10 --lambda 200 --max_iter 130 --fplot out.dat -sigma0 5 -x0 5 -seed 5489
gnuplot -e "filename='out.dat'" cma_multiplt.dem
```
to plot results with matplotlib instead
```
python cma_multiplt.py out.dat
```
to run a check across a range of classical single-objective optimization functions:
```
./test_functions --all
```
for help, do
```
./test_functions --help
```

### Sample code

```C++
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
  std::vector<double> x0(dim,10.0);
  double sigma = 0.1;
  //int lambda = 100; // offsprings at each generation.
  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  //cmaparams.set_algo(BIPOP_CMAES);
  CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
  std::cout << "best solution: " << cmasols << std::endl;
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  return cmasols.run_status();
}
```

### Practical hints

CMA-ES requires two components from the user:
- the initial start point x0;
- the initial value for sigma, the so-called step-size or error guess.

In short: the optimum that is looked after should better not be far away from the interval [x0 - sigma0, x0 + sigma0] in each dimension, where distance is defined by sigma0.

See https://www.lri.fr/~hansen/cmaes_inmatlab.html#practical for more detailed useful advices using CMA-ES.

### Run BBOB 2013 Black-Box Optimization Benchmark

There's an install script in the repository. Do:
```
cd tests
./bbobsetup.sh
```
you can now benchmark any of the implemented flavors of CMA-ES (beware, this make take a while, ~hours):
```
./bbobexperiment -alg bipop
```
for the command above, results will be in repository bipop_bbob
See (7) for more information and details.

### References
- (1) Hansen, N. and A. Ostermeier (2001). Completely Derandomized Self-Adaptation in Evolution Strategies. Evolutionary Computation, 9(2), pp. 159-195.
- (2) Hansen, N. (2009). Benchmarking a BI-Population CMA-ES on the BBOB-2009 Function Testbed. Workshop Proceedings of the GECCO Genetic and Evolutionary Computation Conference, ACM, pp. 2389-2395. http://hal.archives-ouvertes.fr/inria-00382093/en
- (3) N. Hansen, A. Auger, R. Ros, S. Finck, P. Posik: Comparing Results of 31 Algorithms from the Black-Box Optimization Benchmarking BBOB-2009. Workshop Proceedings of the GECCO Genetic and Evolutionary Computation Conference 2010, ACM. http://www.lri.fr/~hansen/ws1p34.pdf
- (4) https://www.lri.fr/~hansen/
- (5) https://www.lri.fr/~hansen/cmaes_inmatlab.html
- (6) Hansen, N., R. Ros (2010). Benchmarking a Weighted Negative Covariance Matrix Update on the BBOB-2010 Noiseless Testbed. Workshop Proceedings of the GECCO Genetic and Evolutionary Computation Conference 2010, ACM, pp. 1673-1680, https://www.lri.fr/~hansen/ws1p32-hansen.pdf
- (7) http://coco.gforge.inria.fr/doku.php?id=bbob-2013
