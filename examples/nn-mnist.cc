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
#include <random>
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
int load_mnist_dataset(const std::string &filename,
		       const int &n,
		       const double &testp,
		       dMat &features, dMat &labels,
		       dMat &tfeatures, dMat &tlabels)
{
  // matrices for features and labels, examples in col
  int trn = ceil(((100.0-testp)/100.0)*n);
  int ttn = n-trn;
  features.resize(784,trn);
  labels.resize(10,trn);
  tfeatures.resize(784,ttn);
  tlabels.resize(10,ttn);
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
      for (size_t i=0;i<strvalues.size();i++)
	values.push_back(strtod(strvalues.at(i).c_str(),NULL));
      if (ne < trn)
	{
	  labels(values.at(0),ne) = 1.0;
	  for (size_t i=1;i<values.size();i++)
	    features(i-1,ne) = values.at(i);
	}
      else
	{
	  tlabels(values.at(0),ne-trn) = 1.0;
	  for (size_t i=1;i<values.size();i++)
	    tfeatures(i-1,ne-trn) = values.at(i);
	}
      ++ne;
    }
  fin.close();
  return 0;
}

// global nn variables etc...
std::vector<int> glsizes;
dMat gfeatures = dMat::Zero(784,100);
dMat glabels = dMat::Zero(10,100);
dMat gtfeatures = dMat::Zero(784,100);
dMat gtlabels = dMat::Zero(10,100);
nn gmnistnn;
bool gsigmoid = false;

// testing
void testing(const CMASolutions &cmasols,
	     const bool &training=true)
{
  dMat cmat = dMat::Zero(10,10);
  Candidate bcand = cmasols.best_candidate();
  std::copy(bcand._x.data(),bcand._x.data()+bcand._x.size(),std::back_inserter(gmnistnn._allparams));
  if (training)
    gmnistnn.forward_pass(gfeatures,glabels);
  else gmnistnn.forward_pass(gtfeatures,gtlabels);
  //std::cout << gmnistnn._lfeatures.cols() << " / " << gmnistnn._lfeatures.rows() << std::endl;
  for (int i=0;i<gmnistnn._lfeatures.cols();i++)
    {
      dMat::Index maxv[2];
      gmnistnn._lfeatures.col(i).maxCoeff(&maxv[0]);
      if (training)
	glabels.col(i).maxCoeff(&maxv[1]);
      else gtlabels.col(i).maxCoeff(&maxv[1]);
      cmat(maxv[1],maxv[0])++;
    }
  std::cerr << "cmat:" << std::endl << cmat << std::endl;
  dMat diago = cmat.diagonal();
  dMat col_sums = cmat.colwise().sum();
  dMat row_sums = cmat.rowwise().sum();
  double precision = diago.transpose().cwiseQuotient(col_sums).sum() / 10.0;
  double recall = diago.cwiseQuotient(row_sums).sum() / 10.0;
  double accuracy = diago.sum() / cmat.sum();
  double f1 = (2 * precision * recall) / (precision + recall);

  if (training)
    std::cout << "training set:\n";
  else std::cout << "testing set:\n";
  std::cout << "precision=" << precision << " / recall=" << recall << std::endl;
  std::cout << "accuracy=" << accuracy << std::endl;
  std::cout << "f1=" << f1 << std::endl;
}

std::random_device rd;
std::mt19937 ggen(rd());
std::uniform_int_distribution<> gunif(0,41999);
int gbatches = -1;

// objective function
FitFunc nn_of = [](const double *x, const int N)
{
  nn mgn = nn(glsizes,gsigmoid);
  for (int i=0;i<N;i++)
    mgn._allparams.push_back(x[i]);
  if (gbatches <= 0)
    mgn.forward_pass(gfeatures,glabels);
  else
    {
      dMat lgfeatures(gfeatures.rows(),gbatches);
      dMat lglabels(glabels.rows(),gbatches);
      for (int j=0;j<gbatches;j++)
	{
	  double u = gunif(ggen);
	  lgfeatures.col(j) = gfeatures.col(u);
	  lglabels.col(j) = glabels.col(u);
	}
      mgn.forward_pass(lgfeatures,lglabels);
    }
  
  //debug
  //std::cout << "net loss= " << gmnistnn._loss << std::endl;
  //debug
  
  return mgn._loss;
};

// gradient function
GradFunc gnn = [](const double *x, const int N)
{
  dVec grad = dVec::Zero(N);
  if (gmnistnn._has_grad)
    {
      gmnistnn._allparams.clear();
      gmnistnn.clear_grad();
      for (int i=0;i<N;i++)
	gmnistnn._allparams.push_back(x[i]);
      gmnistnn.forward_pass(gfeatures,glabels);
      gmnistnn.back_propagate(gfeatures);
      grad = gmnistnn.grad_to_vec(gfeatures.cols());
    }
  //std::cerr << "grad=" << grad.transpose() << std::endl;
  return grad;
};

DEFINE_string(fdata,"train.csv","name of the file that contains the training data for MNIST");
DEFINE_int32(n,100,"max number of examples to train from");
DEFINE_int32(maxsolveiter,-1,"max number of optimization iterations");
DEFINE_string(fplot,"","output file for optimization log");
DEFINE_bool(check_grad,false,"checks on gradient correctness via back propagation");
DEFINE_bool(with_gradient,false,"whether to use the gradient (backpropagation) along with black-box optimization");
DEFINE_double(lambda,-1,"number of offsprings at each generation");
DEFINE_double(sigma0,0.01,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_int32(hlayer,100,"number of neurons in the hidden layer");
DEFINE_bool(sigmoid,false,"whether to use sigmoid units (default is tanh)");
DEFINE_double(testp,0.0,"percentage of the training set used for testing");
DEFINE_int32(mbatch,-1,"size of minibatch");
DEFINE_int32(seed,0,"seed for es");

//TODO: train with batches.
int main(int argc, char *argv[])
{
  ggen.seed(static_cast<uint64_t>(time(nullptr)));
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_check_grad)
    {
      FLAGS_n = 10;
      FLAGS_hlayer = 10;
    }
  if (FLAGS_mbatch > 0)
    gbatches = FLAGS_mbatch;
  
  int err = load_mnist_dataset(FLAGS_fdata,FLAGS_n,FLAGS_testp,gfeatures,glabels,gtfeatures,gtlabels);
  if (err)
    {
      std::cout << "error loading dataset " << FLAGS_fdata << std::endl;
      exit(1);
    }
  gunif = std::uniform_int_distribution<>(0,gfeatures.cols()-1);
  if (FLAGS_check_grad)
    {
      // we check on random features, but we keep the original labels.
      gfeatures.resize(784,FLAGS_n);
      gfeatures = dMat::Random(784,FLAGS_n);
    }
  gsigmoid = FLAGS_sigmoid;
    
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  glsizes = {784, FLAGS_hlayer, 10};
  gmnistnn = nn(glsizes,FLAGS_sigmoid,FLAGS_check_grad || FLAGS_with_gradient);

  if (FLAGS_check_grad)
    {
      if (gmnistnn.grad_check(gfeatures,glabels))
	std::cout << "Gradient check: OK\n";
      else std::cout << "Gradient check did fail\n";
      exit(1);
    }
  
  // training.
  gmnistnn.to_array();
  CMAParameters<> cmaparams(gmnistnn._allparams_dim,&gmnistnn._allparams.front(),FLAGS_sigma0,FLAGS_lambda,FLAGS_seed);
  cmaparams.set_max_iter(FLAGS_maxsolveiter);
  cmaparams._fplot = FLAGS_fplot;
  cmaparams._algo = sepaCMAES;
  cmaparams.set_ftarget(1e-2);
  cmaparams._mt_feval = true;
  /*if (gbatches > 0)
    cmaparams.set_noisy();*/
  CMASolutions cmasols;
  if (!FLAGS_with_gradient)
    cmasols = cmaes<>(nn_of,cmaparams);
  else cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,gnn);
  std::cout << "status: " << cmasols._run_status << std::endl;

  // testing on training set.
  testing(cmasols,true);

  // testing on testing set, if any.
  if (FLAGS_testp)
    testing(cmasols,false);
}
