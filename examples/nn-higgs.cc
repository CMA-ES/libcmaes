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
		       dMat &features, dMat &labels, dMat &weights)
{
  // matrices for features and labels, examples in col 
  features.resize(30,n);
  labels.resize(2,n);
  weights.resize(1,n);
  labels.setZero();
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
      for (size_t i=1;i<values.size()-1;i++) // skipping the event ID.
	features(i-1,ne) = values.at(i);
      std::string label = strvalues.back();
      if (label == "b")
	labels(0,ne) = 1.0;
      else labels(1,ne) = 1.0;
      weights(0,ne) = values.back();
      ++ne;
    }
  fin.close();
  if (n > ne)
    {
      features.resize(30,ne);
      labels.resize(2,ne);
    }
  return 0;
}

// global nn variables etc...
std::vector<int> glsizes;
dMat gfeatures = dMat::Zero(30,100);
dMat glabels = dMat::Zero(2,100);
dMat gweights = dMat::Zero(1,100);
nn ghiggsnn;

double max_ams(const int &n)
{
  double br = 10.0;
  double s = 0.0;
  for (int i=0;i<n;i++)
    {
      if (glabels(1,i) == 1.0) // s
	s += gweights(0,i);
    }
  double cams = sqrt(2.0*((s+br)*log(1.0+s/br)-s));
  //std::cout << "b=" << b << " / s=" << s << " / ams=" << cams << std::endl;
  return -cams;
}

double ams(const nn &hgn)
{
  double br = 10.0;
  double b = 0.0;
  double s = 0.0;
  for (int i=0;i<hgn._lfeatures.cols();i++)
    {
      dMat::Index ind;
      hgn._lfeatures.col(i).maxCoeff(&ind);
      if (ind == 1)
	{
	  if (glabels(ind,i) == 1.0)
	    s += gweights(0,i);
	  else b += gweights(0,i);
	}
    }
  double cams = sqrt(2.0*((s+b+br)*log(1.0+s/(b+br))-s));
  //std::cout << "b=" << b << " / s=" << s << " / ams=" << cams << std::endl;
  return -cams;
}

// objective function
FitFunc nn_of = [](const double *x, const int N)
{
  //std::copy(x,x+N,ghiggsnn._allparams.begin()); // beware.
  /*ghiggsnn._allparams.clear();
  if (ghiggsnn._has_grad)
    ghiggsnn.clear_grad();
  for (int i=0;i<N;i++)
    ghiggsnn._allparams.push_back(x[i]);
    ghiggsnn.forward_pass(gfeatures,glabels);*/

  nn hgn = nn(glsizes);
  for (int i=0;i<N;i++)
    hgn._allparams.push_back(x[i]);
  hgn.forward_pass(gfeatures,glabels);
  
  //debug
  //std::cout << "net loss= " << ghiggsnn._loss << std::endl;
  //debug

  double cams = ams(hgn);
  
  //return ghiggsnn._loss;
  return cams;
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

DEFINE_string(fdata,"train.csv","name of the file that contains the training data for HIGGS");
DEFINE_int32(n,100,"max number of examples to train from");
DEFINE_int32(maxsolveiter,-1,"max number of optimization iterations");
DEFINE_string(fplot,"","output file for optimization log");
DEFINE_bool(check_grad,false,"checks on gradient correctness via back propagation");
DEFINE_bool(with_gradient,false,"whether to use the gradient (backpropagation) along with black-box optimization");
DEFINE_double(lambda,-1,"number of offsprings at each generation");
DEFINE_double(sigma0,1.0,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_int32(hlayer,100,"number of neurons in the hidden layer");

//TODO: train with batches.
int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_check_grad)
    {
      FLAGS_n = 10;
      FLAGS_hlayer = 10;
    }

  int err = load_higgs_dataset(FLAGS_fdata,FLAGS_n,gfeatures,glabels,gweights);
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
    
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  //double minloss = 1e-3; //TODO.
  //int lambda = 1e3;
  glsizes = {30, FLAGS_hlayer, 2};
  ghiggsnn = nn(glsizes,FLAGS_check_grad || FLAGS_with_gradient);

  if (FLAGS_check_grad)
    {
      if (ghiggsnn.grad_check(gfeatures,glabels))
	std::cout << "Gradient check: OK\n";
      else std::cout << "Gradient check did fail\n";
      exit(1);
    }

  std::cout << "max ams=" << -max_ams(FLAGS_n) << std::endl;
  
  // training.
  //double sigma = 2.0;
  std::vector<double> x0(ghiggsnn._allparams_dim,-std::numeric_limits<double>::max());
  CMAParameters<> cmaparams(ghiggsnn._allparams_dim,&x0.front(),FLAGS_sigma0,FLAGS_lambda);
  cmaparams.set_max_iter(FLAGS_maxsolveiter);
  cmaparams._fplot = FLAGS_fplot;
  cmaparams._algo = sepaCMAES;
  //cmaparams.set_ftarget(1e-8);
  CMASolutions cmasols;
  if (!FLAGS_with_gradient)
    cmasols = cmaes<>(nn_of,cmaparams);
  else cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,gnn);
  std::cout << "status: " << cmasols._run_status << std::endl;

  // testing.
  dMat cmat = dMat::Zero(2,2);
  for (int i=0;i<cmasols.best_candidate()._x.size();i++)
    ghiggsnn._allparams.push_back(cmasols.best_candidate()._x(i));
  ghiggsnn.forward_pass(gfeatures,glabels);
  ghiggsnn.forward_pass(gfeatures,glabels);
  //std::cout << ghiggsnn._lfeatures.cols() << " / " << ghiggsnn._lfeatures.rows() << std::endl;
  for (int i=0;i<ghiggsnn._lfeatures.cols();i++)
    {
      dMat::Index maxv[2];
      ghiggsnn._lfeatures.col(i).maxCoeff(&maxv[0]);
      glabels.col(i).maxCoeff(&maxv[1]);
      //std::cout << "maxv=" << maxv << std::endl;
      //if (glabels(maxv[0],i) == 1.0)
      cmat(maxv[1],maxv[0])++;
    }
  //std::cerr << "cmat:" << std::endl << cmat << std::endl;
  dMat diago = cmat.diagonal();
  dMat col_sums = cmat.colwise().sum();
  dMat row_sums = cmat.rowwise().sum();

  /*std::cerr << "col_sums:" << col_sums << std::endl;
    std::cerr << "row_sums:" << row_sums << std::endl;*/
  
  /*double precision = diago.transpose().cwiseQuotient(col_sums).sum() / 10.0;
    double recall = diago.cwiseQuotient(row_sums).sum() / 10.0;*/

  double accuracy = diago.sum() / cmat.sum();
  //double f1 = (2 * precision * recall) / (precision + recall);

  //std::cout << "precision=" << precision << " / recall=" << recall << std::endl;
  std::cout << "accuracy=" << accuracy << std::endl;
  //std::cout << "f1=" << f1 << std::endl;
  std::cout << "max ams=" << -max_ams(FLAGS_n) << std::endl;
  std::cout << "ams=" << -cmasols.best_candidate()._fvalue << std::endl;
}
