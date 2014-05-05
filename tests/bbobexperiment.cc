/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 INRIA
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 * Uses the structure from BBOB.v13.09
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

/* runs an entire experiment for benchmarking MY_OPTIMIZER
* on the noise-free testbed
* or the noisy testbed (change the ifun loop in this case as given below).
*/

#include "cmaes.h"
#include <stdio.h>
#include <string.h>
#include <ctime>
#include <stdlib.h>
#include "bbobStructures.h" /* Include all declarations for BBOB calls */
#include <vector>
#include <gflags/gflags.h>

using namespace libcmaes;

void MY_OPTIMIZER(double(*fitnessfunction)(double*), unsigned int dim, double ftarget, double maxfunevals, int alg)
{
  // map fct to libcmaes FitFunc.
  FitFunc ff = [&](const double *x, const int N)
    {
      return (*fitnessfunction)(const_cast<double*>(x));
    };

  // call to cmaes().
  CMAParameters cmaparams(dim,-1,maxfunevals);
  cmaparams._algo = alg;
  cmaparams._quiet = true;
  cmaes(ff,cmaparams);
  /*CMASolutions cmasols = cmaes(ff,cmaparams);
    Candidate bc = cmasols.best_candidate();*/
}

DEFINE_string(alg,"cmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop & abipop");

int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);

  std::map<int,std::string> flavors;
  if (FLAGS_alg == "cmaes")
    flavors.insert(std::pair<int,std::string>(CMAES_DEFAULT,FLAGS_alg));
  else if (FLAGS_alg == "ipop")
    flavors.insert(std::pair<int,std::string>(IPOP_CMAES,FLAGS_alg));
  else if (FLAGS_alg == "bipop")
    flavors.insert(std::pair<int,std::string>(BIPOP_CMAES,FLAGS_alg));
  else if (FLAGS_alg == "acmaes")
    flavors.insert(std::pair<int,std::string>(aCMAES,FLAGS_alg));
  else if (FLAGS_alg == "aipop")
    flavors.insert(std::pair<int,std::string>(aIPOP_CMAES,FLAGS_alg));
  else if (FLAGS_alg == "abipop")
    flavors.insert(std::pair<int,std::string>(aBIPOP_CMAES,FLAGS_alg));
  
  for (auto mit=flavors.begin();mit!=flavors.end();++mit)
    {
      unsigned int dim[6] = {2, 3, 5, 10, 20, 40};
      unsigned int instances[15] = {1, 2, 3, 4, 5, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
      unsigned int idx_dim, ifun, idx_instances, seed;
      int independent_restarts;
      double maxfunevals, minfunevals;
      
      clock_t t0 = clock();
      time_t Tval;
      /**************************************************
       *          BBOB Mandatory initialization         *
       *************************************************/
      /* retrieve all default parameters of BBOB calls  */
      ParamStruct params = fgeneric_getDefaultPARAMS();
      
      /* modify the following parameters, choosing a different setting
       * for each new experiment */
      std::string reponame = (*mit).second + "_bbob";
      strcpy(params.dataPath,reponame.c_str());  /* different folder for each experiment! */
      /* please beforehand run from the command-line 'python createfolders.py PUT_MY_BBOB_DATA_PATH'
       * to create the necessary folder structure to run an experiment. */
      strcpy(params.algName,(*mit).second.c_str());
      //strcpy(params.comments, "testrun");
      
      /*seed = time(nullptr);
	srand(seed);*/ /* used by MY_OPTIMIZER */
      //printf("random seed set to %d\n", seed);
      
      /* To make the noise deterministic. */
      /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */
      
      /* now the main loop */
      for (idx_dim = 0; idx_dim < 6; idx_dim++)
	{
	  /* Function indices are from 1 to 24 (noiseless) or from 101 to 130 (noisy) */
	  /* for the noisy functions exchange the for loop with */
	  /* for (ifun = 101; ifun <= 130; ifun++) */
	  for (ifun = 1; ifun <= 24; ifun++)
	    {
	      for (idx_instances = 0; idx_instances < 15; idx_instances++)
		{
		  /* set DIM, funcId, instanceId to initialize BBOB fgeneric */
		  params.DIM = dim[idx_dim];
		  params.funcId = ifun;
		  params.instanceId = instances[idx_instances];
		  /* call the BBOB initialization */
		  fgeneric_initialize(params);
		  
		  /* now call your optimizer so that it optimizes the function
		   * fgeneric_evaluate or
		   * fgeneric_evaluate_vector(double * XX, unsigned int howMany,
		   *                          double * result)
		   */
		  
		  /* The fgeneric interface can give some information:
		   *    e.g. fgeneric_ftarget() the target value only for termination
		   *         fgeneric_evaluations() the number of calls to fgeneric_evaluate 
		   *                 after fgeneric_initialization
		   *         fgeneric_best() the best value reached  
		   */
		  maxfunevals = 1e6;//5. * dim[idx_dim]; /* PUT APPROPRIATE MAX. NUMBER OF FEVALS */
		  /* 5. * dim should be fine to just check everything */
		  minfunevals = dim[idx_dim] + 2;  /* PUT MINIMAL USEFUL NUMBER OF FEVALS */
		  independent_restarts = -1;
		  while (fgeneric_evaluations() + minfunevals <= maxfunevals)
		    {
		      if (++independent_restarts > 0) 
                        fgeneric_restart("independent restart");  /* additional info */
		      MY_OPTIMIZER(&fgeneric_evaluate, dim[idx_dim], fgeneric_ftarget(),
				   maxfunevals - fgeneric_evaluations(), (*mit).first);
		      if (fgeneric_best() < fgeneric_ftarget())
                        break;
		    }
		  
		  printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", ifun, dim[idx_dim],
			 instances[idx_instances], fgeneric_evaluations(), independent_restarts);
		  printf(" fbest-ftarget=%.4e, elapsed time [h]: %.2f\n", 
			 fgeneric_best() - fgeneric_ftarget(), (double)(clock()-t0)/CLOCKS_PER_SEC/60./60.);
		  /* call the BBOB closing function to wrap things up neatly */
		  fgeneric_finalize();
		}
	      Tval = time(nullptr);
	      printf("    date and time: %s", ctime(&Tval));
	    }
	  printf("---- dimension %d-D done ----\n", dim[idx_dim]);
	}
    } // end for alg.
    return 0;
}