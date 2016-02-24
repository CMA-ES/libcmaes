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

#include <random>
#include <math.h>

using namespace libcmaes;

void worldFunc(std::vector<double>& world, int & x,int action,std::uniform_real_distribution<> dis,std::mt19937& gen,std::vector<double>& inputs)
{
	//advance in the world if the action is in the direction of the slope
	if ((world[x+1] - world[x] > 0.0) && (action == 1))
	{
		//move one step
		x++;
		//get randomly the new element
		world.push_back(dis(gen));
		//update the input
		inputs[0] = world[x]-world[x-1];
		inputs[1] = world[x+1]-world[x];
	}
	else if ((world[x+1] - world[x] <= 0.0) && (action == -1))
	{
		//move one step
		x++;
		//get randomly the new element
		world.push_back(dis(gen));
		//update the input
		inputs[0] = world[x]-world[x-1];
		inputs[1] = world[x+1]-world[x];
	}
}

void controller(std::vector<double> params, std::vector<double> inputs, int& action)
{
	int idParam = 0;

	double output = 0.0;
	for(unsigned int i = 0 ; i < inputs.size() ; i++)
	{
		output += params[idParam] * inputs[i];
		idParam++;
	}
	output = tanh(output);

	if(output > 0.0)
	{
		action = 1;
	}
	else
	{
		action = -1;
	}

}

FitFunc computeFitness = [](const double *x, const int N)
{
	//sum up the fitnesses computed
  double val = 0.0;
  for (int i=0;i<N;i++)
	{
    val += x[i];
	}

	//turn the problem as a minimization
	//val = 1 / val;
	val = -log(val+1) + log(100);
	std::cout << val << std::endl;
  return val;
};


class customCMAStrategy : public CMAStrategy<CovarianceUpdate>
{
public:
  customCMAStrategy(FitFunc &func,
		    CMAParameters<> &parameters)
    :CMAStrategy<CovarianceUpdate>(func,parameters)
	{
		candidates_uh_set = false;
		nfcalls = 0;
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

		nfcalls = candidates.cols();
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

	bool isset_candidates_uh()
	{
		return candidates_uh_set;
	}
  
  void tell()
  {
    CMAStrategy<CovarianceUpdate>::tell();
  }

  bool stop()
  {
    return CMAStrategy<CovarianceUpdate>::stop();
  }

protected:
	dMat candidates_uh;
  bool candidates_uh_set;
	int nfcalls;
  
};

void evoStep(ESOptimizer<customCMAStrategy,CMAParameters<>>& optim, dMat& candidates, std::vector<RankedCandidate>& nvcandidates, dMat& phenotypes, dMat& phenotypes_uh, int evaluationTime, int popSize, std::vector<double>& params, int& individual, int prevPosition, int position, int timeStep)
{

  if(!optim.stop())
	{
		double fitness = position - prevPosition;

		if(individual < popSize)
		{
			phenotypes((timeStep%evaluationTime),individual) = fitness;
		}
		else
		{
			phenotypes_uh((timeStep%evaluationTime),individual-popSize) = fitness;
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

		if (individual < popSize)
		{
			params.clear();
			for(int i = 0 ; i < candidates.rows() ; i++)
			{
				params.push_back(candidates(i,individual));
			}
		}
		else
		{
			params.clear();
			for(int i = 0 ; i < candidates.rows() ; i++)
			{
				params.push_back(candidates(i,individual-popSize));
			}
		}
	}

	if (optim.stop() == true)
	{
		params.clear();
		dVec tmp = optim.get_solutions().best_candidate().get_x_dvec();
		for(int i = 0 ; i < tmp.cols() ; i++)
		{
			params.push_back(tmp(i));
		}
	}
}

int main(int argc, char *argv[])
{
	//random generator
  std::random_device rd;	
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	//world parameters
	int nbSimulationsStep = 10000000;
	std::vector<double> params ; 
	int prevPosition = 1;
	int position = 1;
	int action = 0;
	std::vector<double> world;
	//generate the elements of the world
	world.push_back(dis(gen));
	world.push_back(dis(gen));
	world.push_back(dis(gen));
	//initialize the inputs
	std::vector<double> inputs;
	inputs.push_back(world[1]-world[0]);
	inputs.push_back(world[2]-world[1]);

	//evolution parameters
  int dim = 2; // problem dimensions.
  std::vector<double> x0(dim,0.5);
  double sigma = 1.0;
	int individual = 0;

	//init evolutionary engine
	std::vector<RankedCandidate> nvcandidates;
  CMAParameters<> cmaparams(x0,sigma);
	cmaparams.set_uh(true);
  ESOptimizer<customCMAStrategy,CMAParameters<>> optim(computeFitness,cmaparams);

	//get the parameters of the evo engine
	int evaluationTime = 100;
	dMat candidates = optim.ask();
	int popSize = candidates.cols();
	dMat phenotypes = dMat(evaluationTime,popSize);
	dMat phenotypes_uh = dMat(evaluationTime,popSize);

	//init the parameters
	params.clear();
	for(int i = 0 ; i < candidates.rows() ; i++)
	{
		params.push_back(candidates(i,individual));
	}

	//run the simulation
	for(int i = 0 ; i < nbSimulationsStep ; i++)
	{
		//get the next move of the agent
		controller(params,inputs,action);

		//update the position of the agent according to its actions
		prevPosition = position;
		worldFunc(world,position,action,dis,gen,inputs);

		//call cma to get the new parameters
		evoStep(optim,candidates,nvcandidates,phenotypes,phenotypes_uh,evaluationTime,popSize,params,individual,prevPosition,position,i);

		if(optim.stop() == true)
		{
			dVec tmp = optim.get_solutions().best_candidate().get_x_dvec();
			for(int i = 0 ; i < tmp.cols() ; i++)
			{
				std::cout << tmp(i) << " " ;
			}
			std::cout << std::endl;
			break;
		}
	}
}

