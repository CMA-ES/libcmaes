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
#include "rankingsvm.hpp"
#include <iostream>

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
  class RSVMSurrogateStrategy : public SimpleSurrogateStrategy<TCovarianceUpdate,TGenoPheno>
  {
  public:
    RSVMSurrogateStrategy()
      :SimpleSurrogateStrategy<TCovarianceUpdate,TGenoPheno>()
    {
    }

    RSVMSurrogateStrategy(FitFunc &func,
			  CMAParameters<TGenoPheno> &parameters)
      :SimpleSurrogateStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters)
    {
      this->_l = 100;
      this->_train = [this](const std::vector<Candidate> &c, const dMat &cov)
	{
	  if (c.empty())
	    return 0;
	  dMat x;
	  dVec fvalues;
	  std::vector<Candidate> cp = c;
	  to_mat_vec(cp,x,fvalues);
	  int niter = 1e6;//floor(50000*sqrt(c.at(0).get_x_size()));
	  dVec xmean = eostrat<TGenoPheno>::get_solutions().xmean();
	  _rsvm = RankingSVM<RBFKernel>();
	  _rsvm._encode = true;
	  _rsvm.train(x,niter,cov,xmean);
	  
	  this->set_train_error(this->compute_error(cp,
						    eostrat<TGenoPheno>::get_solutions().csqinv()));
	  
	  //debug
	  //std::cout << "training error=" << _rsvm.error(x,x,fvalues,cov,xmean) << std::endl;
	  std::cout << "train error=" << this->get_train_error() << std::endl;
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

    ~RSVMSurrogateStrategy() {}
    
    RankingSVM<RBFKernel> _rsvm;
  };

FitFunc fsphere = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
    val += x[i]*x[i];
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
  int dim = 10; // problem dimensions.
  std::vector<double> x0(dim,2.0);
  double sigma = 1.0;

  CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  cmaparams.set_quiet(false);
  cmaparams.set_ftarget(1e-8);
  ESOptimizer<RSVMSurrogateStrategy<>,CMAParameters<>> optim(fsphere,cmaparams);
  optim.set_exploit(false); // test mode.
  
  while(!optim.stop())
    {
      dMat candidates = optim.ask();
      optim.eval(candidates);
      optim.tell();
      optim.inc_iter(); // important step: signals next iteration.
    }
  std::cout << optim.get_solutions() << std::endl;
}
