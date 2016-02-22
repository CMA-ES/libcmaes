/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Author: Jean-Marc Montanier <montanier.jeanmarc@gmail.com>
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
#include <iostream>

using namespace libcmaes;

//TODO do something more complex here that depends on action
void worldFunc(int& x,int action,int& inputs)
{
	x++;
}

//TODO write a little NN implem
void controller(std::vector<double> params, std::vector<double> inputs, int& action)
{
	return 1;
}

FitFunc computeFitness = [](const double *x, const int N)
{
  double val = 0.0;
  for (int i=0;i<N;i++)
	{
    val += x[i];
	}
	val = val / (double)N;
  return val;
};

class customCMAStrategy : public CMAStrategy<CovarianceUpdate>
{
public:
  customCMAStrategy(FitFunc &func,
		    CMAParameters<> &parameters)
    :CMAStrategy<CovarianceUpdate>(func,parameters)
  {
  }

  ~customCMAStrategy() {}

  dMat ask()
  {
    return CMAStrategy<CovarianceUpdate>::ask();
  }

  void eval(const dMat &candidates,
	    const dMat &phenocandidates=dMat(0,0))
  {
    // custom eval.
    for (int r=0;r<candidates.cols();r++)
      {
	_solutions.get_candidate(r).set_x(candidates.col(r));
	if (phenocandidates.size()) // if candidates in phenotype space are given
	  _solutions.get_candidate(r).set_fvalue(_func(phenocandidates.col(r).data(),candidates.rows()));
	else _solutions.get_candidate(r).set_fvalue(_func(candidates.col(r).data(),candidates.rows()));
	
	//std::cerr << "candidate x: " << _solutions.get_candidate(r).get_x_dvec().transpose() << std::endl;
      }

    int nfcalls = candidates.cols();
    // evaluation step of uncertainty handling scheme.
		perform_uh(candidates,phenocandidates,nfcalls);

    update_fevals(nfcalls);
  }
  
  void tell()
  {
    CMAStrategy<CovarianceUpdate>::tell();
  }

  bool stop()
  {
    return CMAStrategy<CovarianceUpdate>::stop();
  }
  
};

int main(int argc, char *argv[])
{
	//world parameters
	int nbSimulationsStep = 1000;
	std::vector<double> params ; 
	int prevPosition = 0;
	int position = 0;
	int action = 0;
	std::vector<double> inputs;
	//TODO initialise the inputs to empty

	//evolution parameters
  int dim = 10; // problem dimensions.
  std::vector<double> x0(dim,10.0);
  double sigma = 0.1;
	int individual = 0;

	//init evolutionary engine
  CMAParameters<> cmaparams(x0,sigma);
  //ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<>> optim(fsphere,cmaparams);
	cmaparams.set_uh(true);
  ESOptimizer<customCMAStrategy,CMAParameters<>> optim(fsphere,cmaparams);

	//TODO init the parameters

	//run the simulation
	for(int i = 0 ; i < nbSimulationsStep ; i++)
	{
		//get the next move of the agent
		controller(params,inputs,action);

		//update the position of the agent according to its actions
		prevPosition = position;
		worldFunc(position,action,inputs);

		//call cma to get the new parameters
		evoStep(optim,params,indivual,prevPosition,position,i);
	}
}

std::vector<double> evoStep(ESOptimizer<customCMAStrategy,CMAParameters<>> optim, std::vector<double>& params, int& individual, int prevPosition, int position, int timeStep)
{
	int evaluationLength = 100;
	dMat candidates = optim.ask();
	int popSize = candidates.cols();
  
  if(!optim.stop())
	{
		double fitness = position - prevPosition;
		//TODO turn this as a minimisation pb

		if(individual < popSize)
		{
			phenotypes(timeStep,individual) = fitness;
		}
		else
		{
			phenotypes_uh(timeStep,individual-popSize) = fitness;
		}

		if((i % evaluationLength) == 0)
		{
			individual++;

			if ((individual == popSize) && (optim.isset_candidates_uh() == false))
			{
				optim.eval(candidates,phenotypes);
				nvcandidates.clear();
			}
		}



		//TODO check that line, it probably does not work. Idea is to get all columns of individual



	}
  std::cout << optim.get_solutions() << std::endl;
}
