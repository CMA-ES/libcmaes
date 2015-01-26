/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
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
#include <mutex>
#include <iostream>

#ifndef GFLAGS_GFLAGS_H_
namespace gflags = google;
#endif  // GFLAGS_GFLAGS_H_

using namespace libcmaes;

void tokenize(const std::string &str,
	      std::vector<std::string> &tokens,
	      const std::string &delim)
{
  
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delim, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delim, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delim, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delim, lastPos);
    }
}

std::mutex fmtx; // WARNING: bbob function calls are NOT thread-safe (learnt the hard way...).

void MY_OPTIMIZER(double(*fitnessfunction)(double*), unsigned int dim, double ftarget, double maxfunevals, int alg, bool noisy, bool withnumgradient, int withtpa)
{
  // map fct to libcmaes FitFunc.
  FitFunc ff = [&](const double *x, const int N)
    {
      std::lock_guard<std::mutex> lck(fmtx);
      double fval = (*fitnessfunction)(const_cast<double*>(x));
      return fval;
    };

  // call to cmaes().
  std::vector<double> x0(dim,-std::numeric_limits<double>::max()); // auto x0 in [-4,4].
  double lbounds[dim];
  double ubounds[dim];
  for (size_t i=0;i<dim;i++)
    {
      lbounds[i] = -5.0;
      ubounds[i] = 5.0;
    }
  GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,dim);
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,2.0,-1,0,gp);
  //CMAParameters<> cmaparams(dim,&x0.front(),2.0,-1,0);
  cmaparams.set_max_fevals(maxfunevals);
  cmaparams.set_ftarget(ftarget);
  //cmaparams.set_x0(-5.0,5.0);
  cmaparams.set_algo(alg);
  cmaparams.set_quiet(true);
  cmaparams.set_gradient(withnumgradient);
  cmaparams.set_tpa(withtpa);
  cmaparams.set_mt_feval(true);
  if (noisy)
    cmaparams.set_noisy();
  cmaes(ff,cmaparams);
  /*CMASolutions cmasols = cmaes(ff,cmaparams);
    Candidate bc = cmasols.best_candidate();
    std::cerr << "solution: " << cmasols << std::endl;*/
}

DEFINE_string(alg,"cmaes","comma separated list of algorithms, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, vdcma");
DEFINE_bool(noisy,false,"whether to benchmark noisy functions");
DEFINE_string(comment,"","comment for the experiment. If using multiple algorithms, the comment will apply to all experiments");
DEFINE_double(maxfunevals,1e6,"maximum number of function evaluations");
DEFINE_double(minfunevals,-1,"minimum number of function evaluations, -1 for automatic definition based on dimension");
DEFINE_bool(with_num_gradient,false,"whether to use numerical gradient injection");
DEFINE_int32(tpa,1,"whether to use two-point adapation for step-size update, 0: no, 1: auto, 2: yes");

int main(int argc, char *argv[])
{
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // parse the alg flags in order to capture all requested algorithm flavors.
  std::vector<std::string> algs;
  tokenize(FLAGS_alg,algs,",");
  
  std::map<int,std::string> flavors;
  for (size_t i=0;i<algs.size();i++)
    {
      if (algs.at(i) == "cmaes")
	flavors.insert(std::pair<int,std::string>(CMAES_DEFAULT,algs.at(i)));
      else if (algs.at(i) == "ipop")
	flavors.insert(std::pair<int,std::string>(IPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "bipop")
	flavors.insert(std::pair<int,std::string>(BIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "acmaes")
	flavors.insert(std::pair<int,std::string>(aCMAES,algs.at(i)));
      else if (algs.at(i) == "aipop")
	flavors.insert(std::pair<int,std::string>(aIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "abipop")
	flavors.insert(std::pair<int,std::string>(aBIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "sepcmaes")
	flavors.insert(std::pair<int,std::string>(sepCMAES,algs.at(i)));
      else if (algs.at(i) == "sepipop")
	flavors.insert(std::pair<int,std::string>(sepIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "sepbipop")
	flavors.insert(std::pair<int,std::string>(sepBIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "sepacmaes")
	flavors.insert(std::pair<int,std::string>(sepaCMAES,algs.at(i)));
      else if (algs.at(i) == "sepaipop")
	flavors.insert(std::pair<int,std::string>(sepaIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "sepabipop")
	flavors.insert(std::pair<int,std::string>(sepaBIPOP_CMAES,algs.at(i)));
      else if (algs.at(i) == "vdcma")
	flavors.insert(std::pair<int,std::string>(VD_CMAES,algs.at(i)));
    }
  
  for (auto mit=flavors.begin();mit!=flavors.end();++mit)
    {
      std::cout << "Running BBOB with algorithm " << (*mit).second << std::endl;
      
      unsigned int dim[6] = {2, 3, 5, 10, 20, 40};
      unsigned int instances[15] = {1, 2, 3, 4, 5, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
      unsigned int idx_dim, ifun, idx_instances;// seed;
      int independent_restarts;
      double maxfunevals = FLAGS_maxfunevals, minfunevals = FLAGS_minfunevals;
      
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
      strcpy(params.comments, FLAGS_comment.c_str());
      
      /*seed = time(nullptr);
	srand(seed);*/ /* used by MY_OPTIMIZER */
      //printf("random seed set to %d\n", seed);
      
      /* To make the noise deterministic. */
      /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */
      
      /* now the main loop */
      for (idx_dim = 0; idx_dim < 6; idx_dim++)
	{
	  /* Function indices are from 1 to 24 (noiseless) or from 101 to 130 (noisy) */
	  unsigned int ifunbegin = 1, ifunend = 24;
	  if (FLAGS_noisy)
	    {
	      ifunbegin = 101;
	      ifunend = 130;
	    }
	  for (ifun = ifunbegin; ifun <= ifunend; ifun++)
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
		  //maxfunevals = 1e6;//5. * dim[idx_dim]; /* PUT APPROPRIATE MAX. NUMBER OF FEVALS */
		  /* 5. * dim should be fine to just check everything */
		  if (minfunevals == -1)
		    minfunevals = dim[idx_dim] + 2;  /* PUT MINIMAL USEFUL NUMBER OF FEVALS */
		  independent_restarts = -1;
		  while (fgeneric_evaluations() + minfunevals <= maxfunevals)
		    {
		      if (++independent_restarts > 0) 
                        fgeneric_restart("independent restart");  /* additional info */
		      MY_OPTIMIZER(&fgeneric_evaluate, dim[idx_dim], fgeneric_ftarget(),
				   maxfunevals - fgeneric_evaluations(), (*mit).first, FLAGS_noisy, FLAGS_with_num_gradient, FLAGS_tpa);
		      if (fgeneric_best() < fgeneric_ftarget())
                        break;
		    }
		  
		  printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", ifun, dim[idx_dim],
			 instances[idx_instances], fgeneric_evaluations(), independent_restarts);
		  //std::cout << "\nfbest=" << fgeneric_best() << " / ftarget=" << fgeneric_ftarget() << std::endl;
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
