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
void worldFunc(int& x,int action,std::vector<int>& inputs)
{
	x++;
}

void controller(std::vector<double> params, std::vector<double> inputs, int& action)
{
	//TODO write a little NN implem
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
	    const dMat &phenocandidates)
	{
		// custom eval.
		for (int r=0;r<candidates.cols();r++)
		{
			_solutions.get_candidate(r).set_x(candidates.col(r));
			_solutions.get_candidate(r).set_fvalue(_func(phenocandidates.col(r).data(),candidates.rows()));
		}

		int nfcalls = candidates.cols();
		// evaluation step of uncertainty handling scheme.

		CMAStrategy<CovarianceUpdate>::select_candidates_uh(candidates,dMat(0,0),candidates_uh);
		candidates_uh_set = true;
	}

	void end_eval(std::vector<RankedCandidate> nvcandidates)
	{
		candidates_uh_set = false;
		CMAStrategy<CovarianceUpdate>::set_candidates_uh(nvcandidates);
		update_fevals(nfcalls);
	}

	dMat get_candidates_uh()
	{
		return candidates_uh;
	}

	void increase_nfcalls()
	{
		nfcalls ++;
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

std::vector<double> evoStep(ESOptimizer<customCMAStrategy,CMAParameters<>> optim, std::vector<double>& params, int& individual, int prevPosition, int position, int timeStep)
{
	int evaluationTime = 100;
	dMat candidates = optim.ask();
	int popSize = candidates.cols();
  
	dMat phenotypes = dMat(evaluationTime,popSize);
	dMat phenotypes_uh = dMat(evaluationTime,popSize);

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

		if((timeStep % evaluationTime) == 0)
		{
			individual++;

			if ((individual == popSize) && (optim.isset_candidates_uh() == false))
			{
				optim.eval(candidates,phenotypes);
				nvcandidates.clear();
			}

			if ((individual == popSize + optim.get_candidates_uh().cols()) && (optim.isset_candidates_uh() == true))
			{

				for (int r=0;r<candidates.cols();r++)
				{
					if (r < optim.get_candidates_uh().cols())
					{
						double nfvalue = computeFitness(phenotypes_uh.col(r).data(),phenotypes_uh.rows());
						Candidate tmp(nfvalue,dVec(candidates.col(r)));
						nvcandidates.emplace_back(nfvalue,tmp,r);
						optim.increase_nfcalls();
					}
					else
					{	
						double nfvalue = computeFitness(phenotypes.col(r).data(),phenotypes_uh.rows());
						Candidate tmp(nfvalue,dVec(candidates.col(r)));
						nvcandidates.emplace_back(nfvalue,tmp,r);
					}
				}

				optim.end_eval(nvcandidates);

				optim.tell();
				optim.inc_iter(); 
				if(!optim.stop())
				{
					candidates = optim.ask();
					individual = 0;
				}
			}
		}

		//TODO check these lines, it probably does not work. Idea is to get all columns of individual
		if (individual < popSize)
		{
			parameters = candidates(,individual);
		}
		else
		{
			parameters = candidates(,individual-popSize);
		}
	}

	if (optim.stop() == true)
	{
		parameters = optim.get_solutions().best_candidate().get_x_dvec();
	}
}

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
	cmaparams.set_uh(true);
  ESOptimizer<customCMAStrategy,CMAParameters<>> optim(computeFitness,cmaparams);

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
		evoStep(optim,params,individual,prevPosition,position,i);
	}
}

