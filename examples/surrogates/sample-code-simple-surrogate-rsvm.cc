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

#include "cmaes.h"
#include "surrogatestrategy.h"
#include "surrogates/rankingsvm.hpp"
#include <iostream>

#include <gflags/gflags.h>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

void to_mat_vec(std::vector<Candidate> &cp,
		dMat &x, dVec &fvalues)
{
  std::sort(cp.begin(),cp.end(),
	    [](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() > c2.get_fvalue();}); // descending sort
  x = dMat(cp.at(0).get_x_size(),cp.size());
  fvalues = dVec(cp.size());
  for (int i=0;i<(int)cp.size();i++)
    {
      x.col(i) = cp.at(i).get_x_dvec().transpose();
      fvalues(i) = cp.at(i).get_fvalue();
    }
}

template <class TGenoPheno> using eostrat = ESOStrategy<CMAParameters<TGenoPheno>,CMASolutions,CMAStopCriteria<TGenoPheno> >;

template<class TCovarianceUpdate=CovarianceUpdate,class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class RSVMSimpleSurrogateStrategy : public SimpleSurrogateStrategy<CMAStrategy,TCovarianceUpdate,TGenoPheno>
  {
  public:
    RSVMSimpleSurrogateStrategy(FitFunc &func,
				CMAParameters<TGenoPheno> &parameters)
      :SimpleSurrogateStrategy<CMAStrategy,TCovarianceUpdate,TGenoPheno>(func,parameters)
    {
      this->_train = [this](const std::vector<Candidate> &c, const dMat &cov)
	{
	  if (c.empty())
	    return 0;
	  dMat x;
	  dVec fvalues;
	  std::vector<Candidate> cp = c;
	  to_mat_vec(cp,x,fvalues);
	  dVec xmean = eostrat<TGenoPheno>::get_solutions().xmean();
	  _rsvm = RankingSVM<RBFKernel>();
	  _rsvm._encode = true;
	  _rsvm.train(x,_rsvm_iter,cov,xmean);
	  
	  dMat cinv = eostrat<TGenoPheno>::get_solutions().csqinv();
	  //this->set_train_error(this->compute_error(cp,cinv)); // clang doesn't like this call within a lambda...
	  this->set_train_error(_rsvm.error(x,x,fvalues,cinv,xmean));
						 	  
	  //debug
	  //std::cout << "training error=" << _rsvm.error(x,x,fvalues,cov,xmean) << std::endl;
	  //std::cout << "train error=" << this->get_train_error() << std::endl;
	  //debug
	  
	  return 0;
	};
      this->_predict = [this](std::vector<Candidate> &c, const dMat &cov)
	{
	  dMat x_test(c.at(0).get_x_size(),c.size());
	  for (int i=0;i<(int)c.size();i++)
	    x_test.col(i) = c.at(i).get_x_dvec().transpose();

	  dMat x_train;
	  dVec fvalues;
	  std::vector<Candidate> tset = this->_tset;
	  to_mat_vec(tset,x_train,fvalues);

	  dVec fit;
	  dVec xmean = eostrat<TGenoPheno>::get_solutions().xmean();
	  _rsvm.predict(fit,x_test,x_train,eostrat<TGenoPheno>::get_solutions().csqinv(),xmean);
	  if (fit.size() != 0)
	    for (int i=0;i<(int)c.size();i++)
	      c.at(i).set_fvalue(fit(i));
	  return 0;
	};
    }

    ~RSVMSimpleSurrogateStrategy() {}
    
    RankingSVM<RBFKernel> _rsvm;
    int _rsvm_iter = 1e6; /**< number of iterations for optimizing the ranking SVM */
};

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
  return val;
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

std::map<std::string,FitFunc> mfuncs;

DEFINE_string(fname,"fsphere","name of the function to optimize");
DEFINE_int32(dim,2,"problem dimension");
DEFINE_int32(lambda,-1,"number of offsprings");
DEFINE_int32(max_iter,-1,"maximum number of iteration (-1 for unlimited)");
DEFINE_int32(max_fevals,-1,"maximum budget as number of function evaluations (-1 for unlimited)");
DEFINE_double(sigma0,-1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_double(ftarget,-std::numeric_limits<double>::infinity(),"objective function target when known");
DEFINE_string(fplot,"","file where to store data for later plotting of results and internal states");
DEFINE_double(x0,-std::numeric_limits<double>::max(),"initial value for all components of the mean vector (-DBL_MAX for automated value)");
DEFINE_bool(no_exploit,false,"whether to exploit the surrogate model");
DEFINE_int32(l,-1,"training set size (number of points)");
DEFINE_int32(rsvm_iter,1e6,"number of iterations for optimizing the ranking SVM");
DEFINE_int32(lifel,-1,"surrogate lifelength, -1 for automatic & dynamic determination");

int main(int argc, char *argv[])
{
  mfuncs["fsphere"]=fsphere;
  mfuncs["elli"]=elli;
  mfuncs["rosenbrock"]=rosenbrock;
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<double> x0(FLAGS_dim,FLAGS_x0);
  
  CMAParameters<> cmaparams(x0,FLAGS_sigma0);
  cmaparams.set_quiet(false);
  cmaparams.set_ftarget(FLAGS_ftarget);
  cmaparams.set_str_algo(FLAGS_alg);
  cmaparams.set_fplot(FLAGS_fplot);
  cmaparams.set_max_iter(FLAGS_max_iter);
  cmaparams.set_max_fevals(FLAGS_max_fevals);
  ESOptimizer<RSVMSimpleSurrogateStrategy<>,CMAParameters<>> optim(fsphere,cmaparams);
  if (FLAGS_no_exploit)
    optim.set_exploit(!FLAGS_no_exploit);
  optim.set_nsteps(FLAGS_lifel);
  optim._rsvm_iter = FLAGS_rsvm_iter;
  
  while(!optim.stop())
    {
      dMat candidates = optim.ask();
      optim.eval(candidates);
      optim.tell();
      optim.inc_iter(); // important step: signals next iteration.
    }
  std::cout << optim.get_solutions() << std::endl;
}
