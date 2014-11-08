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

#include "surrogates/rankingsvm.hpp"
#include "cmaes.h"
#include <chrono>

using namespace libcmaes;

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

void sample_points(const int &n,
		   const int &d,
		   dMat &x,
		   dVec &fvalues)
{
  std::vector<Candidate> c;
  for (int i=0;i<n;i++)
    {
      dVec x = dVec::Random(d)*10.0; // in [-10,10]
      double fvalue = rosenbrock(x.data(),d);
      c.push_back(Candidate(fvalue,x));
    }
  std::sort(c.begin(),c.end(),
	    [](Candidate const &c1, Candidate const &c2){return c1.get_fvalue() > c2.get_fvalue();}); // ranking in descending order.
  x = dMat(d,n);
  fvalues = dVec(n);
  for (int i=0;i<n;i++)
    {
      x.col(i) = c.at(i).get_x_dvec().transpose();
      fvalues(i) = c.at(i).get_fvalue();
    }
}

// generate data.
const int n = 1000;
const int nt = 100;
const int d = 10;
const int niter = 100000;

int main()
{
  dMat x;
  dVec fvalues;
  sample_points(n,d,x,fvalues);
  //RankingSVM<PolyKernel<2>> rsvm;
  RankingSVM<RBFKernel> rsvm;
  dMat cov = dMat(0,0);
  dVec xm = dVec(d); // unused xmean.
  std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
  rsvm.train(x,niter,cov,xm);
  std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
  std::cout << "training time=" << std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count() / 1000.0 << " seconds\n";
  std::cout << "training error=" << rsvm.error(x,x,fvalues,cov,xm) << std::endl;
  dMat tx;
  dVec tfvalues;
  sample_points(nt,d,tx,tfvalues);
  std::cout << "test error=" << rsvm.error(tx,x,tfvalues,cov,xm) << std::endl;
}
