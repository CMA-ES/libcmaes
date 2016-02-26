/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Demonstration of online optimisation
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
	if (((world[x+1] - world[x] > 0.0) && (action == 1)) || ((world[x+1] - world[x] <= 0.0) && (action == -1)))
	{
		//move one step
		x++;

		//get randomly the new element
		double propNewPos = dis(gen);
		//generate position until the difference is big enough
		while(fabs(propNewPos - world[x]) < 0.05)
		{
			propNewPos = dis(gen);
		}
		world.push_back(propNewPos);

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
	val = -log(val+1);
  return val;
};


class onlineCMAStrategy : public CMAStrategy<CovarianceUpdate>
{
public:
  onlineCMAStrategy(FitFunc &func,
		    CMAParameters<> &parameters)
    :CMAStrategy<CovarianceUpdate>(func,parameters)
	{
		candidates_uh_set = false;
		nfcalls = 0;
		genotype = 0;

	}

  ~onlineCMAStrategy() {}

	void init(int inEvaluationTime)
	{
		evaluationTime = inEvaluationTime;
		//init the parameters of the evo engine
		candidates = ask();
		popSize = candidates.cols();
		phenotypes = dMat(evaluationTime,popSize);
		phenotypes_uh = dMat(evaluationTime,popSize);

		nvcandidates.clear();
	}

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

	void setFitness(int timeStep, double fitness)
	{
		if(genotype < popSize)
		{
			phenotypes((timeStep%evaluationTime),genotype) = fitness;
		}
		else
		{
			phenotypes_uh((timeStep%evaluationTime),(genotype-popSize)) = fitness;
		}
	}

	void completeEvaluation()
	{
		genotype++;

		//if we are at the end of the regular evaluation, perform the evals
		//eval will set the candidates_uh
		if ((genotype == popSize) && (candidates_uh_set == false))
		{
			eval(candidates,phenotypes);
			nvcandidates.clear();
		}

		//if we are at the end of the re-evaluations we have to  place the result
		//and get the new generation
		if ((genotype == popSize + candidates_uh.cols()) && (candidates_uh_set == true))
		{
			for (int r=0;r<candidates.cols();r++)
			{
				//place the result of  the re-evaluation of re-evaluated genotypes
				if (r < candidates_uh.cols())
				{
					double nfvalue = computeFitness(phenotypes_uh.col(r).data(),phenotypes_uh.rows());
					Candidate tmp(nfvalue,dVec(candidates.col(r)));
					nvcandidates.emplace_back(nfvalue,tmp,r);
					increase_nfcalls();
				}
				//place the results previously obtained for the other genotypes
				else
				{	
					double nfvalue = computeFitness(phenotypes.col(r).data(),phenotypes_uh.rows());
					Candidate tmp(nfvalue,dVec(candidates.col(r)));
					nvcandidates.emplace_back(nfvalue,tmp,r);
				}
			}

			//finish the evaluation
			end_eval(nvcandidates);
			tell();
			inc_iter(); 

			//if we have not finish prepare the next generation
			if(!stop())
			{
				candidates = ask();
				genotype = 0;
			}
		}
	}

	std::vector<double> getParams()
	{
		std::vector<double> params;
		if (stop() == true)
		{
			dVec tmp = get_solutions().best_candidate().get_x_dvec();
			for(int i = 0 ; i < tmp.cols() ; i++)
			{
				params.push_back(tmp(i));
			}
		}
		else
		{
			if (genotype < popSize)
			{
				params.clear();
				for(int i = 0 ; i < candidates.rows() ; i++)
				{
					params.push_back(candidates(i,genotype));
				}
			}
			else
			{
				params.clear();
				for(int i = 0 ; i < candidates.rows() ; i++)
				{
					params.push_back(candidates(i,genotype-popSize));
				}
			}
		}
		return params;
	}

protected:
	dMat candidates;
	dMat candidates_uh;
  bool candidates_uh_set;
	int nfcalls;

	int genotype;
	int popSize;
	int evaluationTime;
	
	dMat phenotypes;
	dMat phenotypes_uh;

  std::vector<RankedCandidate> nvcandidates;
};

void evoStep(ESOptimizer<onlineCMAStrategy,CMAParameters<>>& optim, int evaluationTime, std::vector<double>& params, int prevPosition, int position, int timeStep)
{

  if(!optim.stop())
	{
		//compute the fitness for this timestep
		double fitness = position - prevPosition;
		
		//set the fitness
		optim.setFitness(timeStep, fitness);

		//if we have reach  the end of evaluation time, ask to CMA-ES to do its job
		if((timeStep % evaluationTime) == 0)
		{
			optim.completeEvaluation();
		}
	}

	//get the parameters from CMA-ES
	params = optim.getParams();

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
	//generate the three first elements of the world
	world.push_back(dis(gen));
	world.push_back(dis(gen));
	world.push_back(dis(gen));
	//initialize the inputs of the controller
	std::vector<double> inputs;
	inputs.push_back(world[1]-world[0]);
	inputs.push_back(world[2]-world[1]);

	//evolution parameters
  int dim = 2; // problem dimensions.
  std::vector<double> x0(dim,0.5);
  double sigma = 1.0;
	int evaluationTime = 100;

	//initialize CMA-ES
  CMAParameters<> cmaparams(x0,sigma);
	cmaparams.set_uh(true);
  ESOptimizer<onlineCMAStrategy,CMAParameters<>> optim(computeFitness,cmaparams);
	optim.init(evaluationTime);

	//initalize the parameters of the controller
	params = optim.getParams();

	//run the simulation
	for(int i = 0 ; i < nbSimulationsStep ; i++)
	{
		//get the next move of the agent
		controller(params,inputs,action);

		//update the position of the agent according to its actions
		//this will generate the next world element if needed, i.e. if the agent is moving
		prevPosition = position;
		worldFunc(world,position,action,dis,gen,inputs);

		//Run the evolutionary step
		evoStep(optim,evaluationTime,params,prevPosition,position,i);

		//if the optimization is done, display the result
		if(optim.stop() == true)
		{
			std::cout << "stop" << std::endl;
			dVec tmp = optim.get_solutions().best_candidate().get_x_dvec();
			for(int i = 0 ; i < tmp.rows() ; i++)
			{
				std::cout << tmp(i) << " " ;
			}
			std::cout << std::endl;
			break;
		}
	}
}

