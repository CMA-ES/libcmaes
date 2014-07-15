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

#include "nn.h"
#include "cmaes.h"
#include <fstream>
#include <gflags/gflags.h>
#include <cstdlib>
#include <chrono>
#include <iostream>

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

// load dataset into memory
int load_higgs_dataset(const std::string &filename,
		       const int &n,
		       const int &testp,
		       dMat &features, dMat &labels, dMat &weights,
		       dMat &tfeatures, dMat &tlabels, dMat &tweights)
{
  // matrices for features and labels, examples in col
  int trn = ceil(((100.0-testp)/100.0)*n);
  int ttn = n-trn;
  features.resize(30,trn);
  labels.resize(2,trn);
  weights.resize(1,trn);
  tfeatures.resize(30,ttn);
  tlabels.resize(2,ttn);
  tweights.resize(1,ttn);
  labels.setZero();
  tlabels.setZero();
  std::ifstream fin(filename);
  if (!fin.good())
    return -1;
  std::string line;
  std::getline(fin,line); // bypass first line
  int ne = 0;
  while(std::getline(fin,line)
	&& ne < n)
    {
      //std::cout << "line: " << line << std::endl;
      std::vector<std::string> strvalues;
      tokenize(line,strvalues,",");
      std::vector<double> values;
      for (size_t i=0;i<strvalues.size()-1;i++)
	{
	  double val = strtod(strvalues.at(i).c_str(),NULL);
	  if (val == -999.0)
	    val = 0.0; // XXX: not sure what to do with the out of bounds (by dataset definition) values.
	  values.push_back(val);
	}
      std::string label = strvalues.back();
      if (ne < trn)
	{
	  for (size_t i=1;i<values.size()-1;i++) // skipping the event ID.
	    features(i-1,ne) = values.at(i);
	  if (label == "b")
	    labels(0,ne) = 1.0;
	  else labels(1,ne) = -1.0;
	  weights(0,ne) = values.back();
	}
      else
	{
	  for (size_t i=1;i<values.size()-1;i++) // skipping the event ID.
	    tfeatures(i-1,ne-trn) = values.at(i);
	  if (label == "b")
	    tlabels(0,ne-trn) = 1.0;
	  else tlabels(1,ne-trn) = -1.0;
	  tweights(0,ne-trn) = values.back();
	}
      ++ne;
    }
  fin.close();
  return 0;
}

// global nn variables etc...
std::vector<int> glsizes;
dMat gfeatures = dMat::Zero(30,100);
dMat glabels = dMat::Zero(2,100);
dMat gweights = dMat::Zero(1,100);
dMat gtfeatures = dMat::Zero(30,100);
dMat gtlabels = dMat::Zero(2,100);
dMat gtweights = dMat::Zero(1,100);
nn ghiggsnn;
bool gsigmoid = false;
bool gtraining = true;
bool ggradient = false;

double max_ams(const int &n,
	       const bool &training=true)
{
  double br = 10.0;
  double s = 0.0;
  if (training)
    {
      for (int i=0;i<glabels.cols();i++)
	{
	  if (glabels(1,i) == -1.0) // s
	    s += gweights(0,i);
	}
    }
  else
    {
      for (int i=0;i<gtlabels.cols();i++)
	{
	  if (gtlabels(1,i) == -1.0) // s
	    s += gtweights(0,i);
	}
    }
  double cams = sqrt(2.0*((s+br)*log(1.0+s/br)-s));
  //std::cout << "b=" << b << " / s=" << s << " / ams=" << cams << std::endl;
  return -cams;
}

double ams(const nn &hgn,
	   const dMat &gweights,
	   const dMat &glabels)
{
  double br = 10.0;
  double b = 0.0;
  double s = 0.0;
  //std::cerr << "lfeatures:" << hgn._lfeatures << std::endl;
  for (int i=0;i<hgn._lfeatures.cols();i++)
    {
      dMat::Index ind;
      hgn._lfeatures.col(i).cwiseAbs().maxCoeff(&ind);
      if (ind == 1)
	{
	  if (glabels(ind,i) == -1.0)
	    s += gweights(0,i);
	  else b += gweights(0,i);
	}
    }
  double cams = sqrt(2.0*((s+b+br)*log(1.0+s/(b+br))-s));
  //std::cout << "b=" << b << " / s=" << s << " / ams=" << cams << std::endl;
  return -cams;
}

// testing
void testing(const CMASolutions &cmasols,
	     const bool &training=true)
{
  std::cout << "testing / training=" << training << std::endl;
  gtraining = training;
  dMat cmat = dMat::Zero(2,2);
  ghiggsnn._allparams.clear();
  for (int i=0;i<cmasols.best_candidate()._x.size();i++)
    ghiggsnn._allparams.push_back(cmasols.best_candidate()._x(i));
  if (training)
    ghiggsnn.forward_pass(gfeatures,glabels);
  else ghiggsnn.forward_pass(gtfeatures,gtlabels);
  //std::cout << ghiggsnn._lfeatures.cols() << " / " << ghiggsnn._lfeatures.rows() << std::endl;
  for (int i=0;i<ghiggsnn._lfeatures.cols();i++)
    {
      dMat::Index maxv[2];
      ghiggsnn._lfeatures.col(i).cwiseAbs().maxCoeff(&maxv[0]);
      if (training)
	glabels.col(i).cwiseAbs().maxCoeff(&maxv[1]);
      else gtlabels.col(i).cwiseAbs().maxCoeff(&maxv[1]);
      cmat(maxv[1],maxv[0])++;
    }
  std::cerr << "cmat:" << std::endl << cmat << std::endl;
  dMat diago = cmat.diagonal();
  dMat col_sums = cmat.colwise().sum();
  dMat row_sums = cmat.rowwise().sum();
  double precision = diago.transpose().cwiseQuotient(col_sums).sum() / 2.0;
  double recall = diago.cwiseQuotient(row_sums).sum() / 2.0;
  double accuracy = diago.sum() / cmat.sum();
  double f1 = (2 * precision * recall) / (precision + recall);

  if (training)
    std::cout << "\n**** training set (" << gfeatures.cols() << "):\n";
  else std::cout << "\n**** testing set (" << gtfeatures.cols() << "):\n";
  std::cout << "precision=" << precision << " / recall=" << recall << " / accuracy=" << accuracy << " / f1=" << f1 << std::endl;
  std::cout << "max ams=" << -max_ams(gfeatures.cols(),training) << " / ams=" << -ams(ghiggsnn,training ? gweights : gtweights, training ? glabels : gtlabels) << std::endl;
}

// objective function
FitFunc nn_of = [](const double *x, const int N)
{
  nn hgn = nn(glsizes,gsigmoid);
  for (int i=0;i<N;i++)
    hgn._allparams.push_back(x[i]);
  hgn.forward_pass(gfeatures,glabels);
  
  //debug
  /*if (ggradient)
    std::cout << "net loss= " << hgn._loss << std::endl;*/
  //debug

  if (!ggradient)
    {
      double cams = ams(hgn,gtraining ? gweights : gtweights, gtraining ? glabels : gtlabels);
      return cams;
    }
  else return hgn._loss;
};

// gradient function
GradFunc gnn = [](const double *x, const int N)
{
  dVec grad = dVec::Zero(N);
  if (ghiggsnn._has_grad)
    {
      ghiggsnn._allparams.clear();
      ghiggsnn.clear_grad();
      for (int i=0;i<N;i++)
	ghiggsnn._allparams.push_back(x[i]);
      ghiggsnn.forward_pass(gfeatures,glabels);
      ghiggsnn.back_propagate(gfeatures);
      grad = ghiggsnn.grad_to_vec(gfeatures.cols());
    }
  //std::cerr << "grad=" << grad.transpose() << std::endl;
  return grad;
};

ProgressFunc<CMAParameters<>,CMASolutions> hpfunc = [](const CMAParameters<> &cmaparams, const CMASolutions &cmasols)
{
  nn hgn = nn(glsizes,gsigmoid);
  for (int i=0;i<cmaparams._dim;i++)
    hgn._allparams.push_back(cmasols.best_candidate()._x(i));
  hgn.forward_pass(gfeatures,glabels);
  double cams = ams(hgn,gtraining ? gweights : gtweights, gtraining ? glabels : gtlabels);
  std::cout << "iter=" << cmasols._niter << " / evals=" << cmaparams._lambda * cmasols._niter << " / f-value=" << cmasols._best_candidates_hist.back()._fvalue <<  " / ams=" << cams << " / sigma=" << cmasols._sigma << (cmaparams._lazy_update && cmasols._updated_eigen ? " / cupdate="+std::to_string(cmasols._updated_eigen) : "") << " " << cmasols._elapsed_last_iter << std::endl;
  return 0;
};

DEFINE_string(fdata,"train.csv","name of the file that contains the training data for HIGGS");
DEFINE_int32(n,100,"max number of examples to train from");
DEFINE_int32(maxsolveiter,-1,"max number of optimization iterations");
DEFINE_string(fplot,"","output file for optimization log");
DEFINE_bool(check_grad,false,"checks on gradient correctness via back propagation");
DEFINE_bool(with_gradient,false,"whether to use the gradient (backpropagation) along with black-box optimization");
DEFINE_int32(lambda,-1,"number of offsprings at each generation");
DEFINE_double(sigma0,1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_string(hlayers,"100","comma separated list of number of neurons per hidden layer");
DEFINE_bool(sigmoid,false,"whether to use sigmoid units (default is tanh)");
DEFINE_double(testp,0.0,"percentage of the training set used for testing");
DEFINE_string(alg,"sepacmaes","algorithm, among cmaes, ipop, bipop, acmaes, aipop, abipop, sepcmaes, sepipop, sepbipop, sepacmaes, sepaipop, sepabipop");
DEFINE_double(lbound,std::numeric_limits<double>::max()/-1e2,"lower bound to parameter vector");
DEFINE_double(ubound,std::numeric_limits<double>::max()/1e2,"upper bound to parameter vector");
DEFINE_double(x0,-std::numeric_limits<double>::max(),"initial value for all components of the mean vector (-DBL_MAX for automated value)");
DEFINE_bool(mbatch,false,"whether to use minibatches");
DEFINE_int32(mbatch_budget,-1,"max budget when using minibatches");
DEFINE_double(mbatch_ftarget,3.0,"AMS target when using minibatches");
DEFINE_bool(mbatch_sim,false,"simplified output for minibatches in order to pipe to file");
DEFINE_bool(with_num_gradient,false,"whether to use numerical gradient injection");

//TODO: train with batches.
int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_check_grad)
    {
      FLAGS_n = 10;
      FLAGS_hlayers = "10";
    }

  //FLAGS_n = 10000;
  std::vector<std::string> hlayers_str;
  std::vector<int> hlayers;
  tokenize(FLAGS_hlayers,hlayers_str,",");
  for (size_t i=0;i<hlayers_str.size();i++)
    hlayers.push_back(atoi(hlayers_str.at(i).c_str()));

  int load_size = FLAGS_n;
  if (FLAGS_mbatch)
    load_size = 250000;
  int err = load_higgs_dataset(FLAGS_fdata,load_size,FLAGS_testp,gfeatures,glabels,gweights,gtfeatures,gtlabels,gtweights);
  if (err)
    {
      std::cout << "error loading dataset " << FLAGS_fdata << std::endl;
      exit(1);
    }
  
  if (FLAGS_check_grad)
    {
      // we check on random features, but we keep the original labels.
      gfeatures.resize(30,FLAGS_n);
      gfeatures = dMat::Random(30,FLAGS_n);
    }
  gsigmoid = FLAGS_sigmoid;
  
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  glsizes.push_back(30);
  for (size_t i=0;i<hlayers.size();i++)
    glsizes.push_back(hlayers.at(i));
  glsizes.push_back(2);
  ghiggsnn = nn(glsizes,gsigmoid,FLAGS_check_grad || FLAGS_with_gradient);
  ggradient = FLAGS_with_gradient;
  
  if (FLAGS_check_grad)
    {
      if (ghiggsnn.grad_check(gfeatures,glabels))
	std::cout << "Gradient check: OK\n";
      else std::cout << "Gradient check did fail\n";
      exit(1);
    }

  std::cout << "max ams=" << -max_ams(gfeatures.cols()) << std::endl;
  
  // training.
  /*double lbounds[ghiggsnn._allparams_dim];
  double ubounds[ghiggsnn._allparams_dim];
  for (size_t i=0;i<ghiggsnn._allparams_dim;i++)
    {
      lbounds[i] = FLAGS_lbound;
      ubounds[i] = FLAGS_ubound;
    }
    GenoPheno<pwqBoundStrategy,NoScalingStrategy> gp(lbounds,ubounds,ghiggsnn._allparams_dim);*/

  dMat ggfeatures;
  dMat gglabels;
  dMat ggweights;
  CMASolutions cmasols;
  int npasses = 1;
  if (FLAGS_mbatch)
    npasses = ceil(gfeatures.cols()/static_cast<double>(FLAGS_n));
  std::vector<double> sigma0(npasses,FLAGS_sigma0);
  double gams = 0.0;
  int nevals = 0;
  int elapsed = 0;
  int elapsed_total = 0;
  bool init = false;
  bool run = true;

  std::cout << "npasses=" << npasses << std::endl;
  std::cout << "dim=" << ghiggsnn._allparams_dim << std::endl;

  std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
  while (run)
    {
      for (int i=0;i<npasses;i++)
	{
	  if (gams >= FLAGS_mbatch_ftarget
	      || (FLAGS_mbatch_budget != -1 && nevals >= FLAGS_mbatch_budget))
	    {
	      run = false;
	      break;
	    }
	  std::vector<double> x0;
	  if (i == 0 && !init)
	    {
	      ggfeatures = gfeatures;
	      gglabels = glabels;
	      ggweights = gweights;
	      ghiggsnn.to_array();
	      if (FLAGS_x0==-std::numeric_limits<double>::max())
		x0 = ghiggsnn._allparams;
	      else x0 = std::vector<double>(ghiggsnn._allparams_dim,FLAGS_x0);
	      init = true;
	    }
	  else
	    {
	      
	      Candidate bcand = cmasols.best_candidate();
	      std::copy(bcand._x.data(),bcand._x.data()+bcand._x.size(),std::back_inserter(x0));
	      nn hgn = nn(glsizes,gsigmoid);
	      for (int i=0;i<(int)x0.size();i++)
		hgn._allparams.push_back(x0[i]);
	      hgn.forward_pass(ggfeatures,gglabels);
	      gams = -ams(hgn,ggweights,gglabels);
	    }

	  std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	  elapsed_total = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	  if (!FLAGS_mbatch_sim)
	    std::cout << "pass #" << i << " / ams=" << gams << " / nevals=" << nevals << " / tim=" << elapsed/1000.0 << " / timt=" << elapsed_total/1000.0 << std::endl;
	  else std::cout << gams << "," << nevals << "," << elapsed/1000.0 << "\t" << elapsed_total / 1000.0 << std::endl;
	  
	  int beg = i*FLAGS_n;
	  int bsize = FLAGS_n;
	  if (i == npasses-1)
	    bsize = ggfeatures.cols()-i*FLAGS_n;
	  gfeatures = ggfeatures.block(0,beg,ggfeatures.rows(),bsize);
	  glabels = gglabels.block(0,beg,gglabels.rows(),bsize);
	  gweights = ggweights.block(0,beg,ggweights.rows(),bsize);
	  CMAParameters<> cmaparams(ghiggsnn._allparams_dim,&x0.front(),sigma0[i],FLAGS_lambda);//,0,gp);
	  cmaparams.set_max_iter(FLAGS_maxsolveiter);
	  cmaparams._fplot = FLAGS_fplot;
	  cmaparams._with_gradient = FLAGS_with_num_gradient;
	  //cmaparams._lazy_update = true;
	  if (FLAGS_alg == "cmaes")
	    cmaparams._algo = CMAES_DEFAULT;
	  else if (FLAGS_alg == "ipop")
	    cmaparams._algo = IPOP_CMAES;
	  else if (FLAGS_alg == "bipop")
	    cmaparams._algo = BIPOP_CMAES;
	  else if (FLAGS_alg == "acmaes")
	    cmaparams._algo = aCMAES;
	  else if (FLAGS_alg == "aipop")
	    cmaparams._algo = aIPOP_CMAES;
	  else if (FLAGS_alg == "abipop")
	    cmaparams._algo = aBIPOP_CMAES;
	  else if (FLAGS_alg == "sepcmaes")
	    cmaparams._algo = sepCMAES;
	  else if (FLAGS_alg == "sepipop")
	    cmaparams._algo = sepIPOP_CMAES;
	  else if (FLAGS_alg == "sepbipop")
	    cmaparams._algo = sepBIPOP_CMAES;
	  else if (FLAGS_alg == "sepacmaes")
	    cmaparams._algo = sepaCMAES;
	  else if (FLAGS_alg == "sepaipop")
	    cmaparams._algo = sepaIPOP_CMAES;
	  else if (FLAGS_alg == "sepabipop")
	    cmaparams._algo = sepaBIPOP_CMAES;
	  else
	    {
	      std::cout << "unknown algorithm flavor " << FLAGS_alg << std::endl;
	      exit(-1);
	    }
	  //cmaparams.set_ftarget(1e-8);
	  cmaparams._mt_feval = true;
	  if (FLAGS_mbatch)
	    cmaparams._quiet = true;
	  //if (!FLAGS_with_gradient)
	  cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,nullptr,cmasols);
	  sigma0[i] = cmasols._sigma;
	  nevals += cmasols._nevals;
	  elapsed += cmasols._elapsed_time;
	  //else cmasols = cmaes<>(nn_of,cmaparams,hpfunc,gnn);
	  //std::cout << "status: " << cmasols._run_status << std::endl;
	}
      if (!FLAGS_mbatch)
	break;
    }
  
  gfeatures = ggfeatures;
  glabels = gglabels;
  gweights = ggweights;
  
  // testing on training set
  testing(cmasols,true);
  
  // testing on test set
  if (FLAGS_testp)
    testing(cmasols,false);
}
  
