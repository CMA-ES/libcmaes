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

#include "cmaes.h"
#include <fstream>
#include <glog/logging.h>
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
int load_mnist_dataset(const std::string &filename,
		       const int &n,
		       dMat &features, dMat &labels)
{
  // matrices for features and labels, examples in col 
  features.resize(784,n);
  labels.resize(10,n);
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
      for (size_t i=0;i<strvalues.size();i++)
	values.push_back(strtod(strvalues.at(i).c_str(),NULL));
      labels(values.at(0),ne) = 1.0;
      for (size_t i=1;i<values.size();i++)
	features(i-1,ne) = values.at(i);
      ++ne;
    }
  fin.close();
  if (n > ne)
    {
      features.resize(784,ne);
      labels.resize(1,ne);
    }
  return 0;
}

// build the network
//- input layer
//- hidden layer
//- top layer + loss function.
class nn
{
public:
  nn()
  {};
  
  nn(const std::vector<int> &lsizes)
    :_lsizes(lsizes)
  {
    for (size_t i=0;i<_lsizes.size()-1;i++)
      {
	_lweights.push_back(dMat::Random(_lsizes.at(i),_lsizes.at(i+1)));
	_lb.push_back(dVec::Zero(_lsizes.at(i+1)));
	_allparams_dim += _lweights.at(i).size() + _lsizes.at(i);
      }
  }
  
  ~nn() {};

  static dMat sigmoid(const dMat &M)
  {
    dMat expM(M.rows(),M.cols());
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    expM(i,j) = exp(-M(i,j));
	  }
      }
    dMat denom = expM + dMat::Constant(M.rows(),M.cols(),1);
    return denom.cwiseInverse();
  }

  static dMat softmax(const dMat &M)
  {
    dMat expM(M.rows(),M.cols());
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    expM(i,j) = exp(M(i,j));
	  }
      }
    dMat sums = expM.colwise().sum();
    // div by row vector.
    for (int i=0;i<expM.cols();i++)
      {
	for (int j=0;j<expM.rows();j++)
	  {
	    expM(j,i) /= sums(0,i);
	  }
      }
    return expM;
  }

  dMat get_loss(const dMat &prediction, const dMat &labels)
  {
    dMat M = labels.cwiseProduct(prediction);
    dMat Mc = M.colwise().sum();

    dMat logM(Mc.rows(),Mc.cols());
    for (int i=0;i<Mc.rows();i++)
      {
	for (int j=0;j<Mc.cols();j++)
	  {
	    logM(i,j) = Mc(i,j) == 0 ? 0 : log(Mc(i,j));
	  }
      }
    return (-1) * logM;
  }
    
  // forward pass over all features at once.
  void forward_pass(const dMat &features,
		    const dMat &labels)
  {
    to_matrices();
    for (size_t i=0;i<_lweights.size();i++)
      {
	dMat activation;
	if (i == 0)
	  activation = (_lweights.at(i).transpose() * features).colwise() + _lb.at(i);
	else activation = (_lweights.at(i).transpose() * _lfeatures).colwise() + _lb.at(i);
	if (i != _lweights.size()-1)
	  _lfeatures = sigmoid(activation);
	else _lfeatures = softmax(activation);
      }
    
    // loss.
    if (labels.size() > 0) // training mode.
      {
	//dMat delta = _lfeatures - labels;
	/*std::cout << "features:\n";
	std::cout << lfeatures << std::endl;
	std::cout << "labels:\n";
	std::cout << labels << std::endl;
	std::cout << "delta:\n";
	std::cout << delta << std::endl;*/
	//_loss = delta.norm();
	_loss = get_loss(_lfeatures,labels).mean();
	//std::cerr << "loss=" << _loss << std::endl;
      }
  };

  void to_array()
  {
    _allparams.clear();
    _allparams.reserve(_allparams_dim);
    auto vit = _allparams.begin();
    for (size_t i=0;i<_lweights.size();i++)
      {
	std::copy(_lweights.at(i).data(),_lweights.at(i).data()+_lweights.at(i).size(),vit);
	vit += _lweights.at(i).size();
	std::copy(_lb.at(i).data(),_lb.at(i).data()+_lsizes.at(i+1),vit);
	vit += _lsizes.at(i+1);
      }
  }

  void to_matrices()
  {
    auto vit = _allparams.begin();
    for (size_t i=0;i<_lweights.size();i++)
      {
	std::copy(vit,vit+_lweights.at(i).size(),_lweights.at(i).data());
        vit += _lweights.at(i).size();
	std::copy(vit,vit+_lsizes.at(i+1),_lb.at(i).data());
	vit += _lsizes.at(i+1);
      }
  }
  
  std::vector<int> _lsizes; /**< layer sizes. */
  std::vector<dMat> _lweights; /**< weight matrice, per layer. */
  std::vector<dVec> _lb;
  unsigned int _allparams_dim = 0;
  std::vector<double> _allparams; /**< all parameters, flat representation. */
  dMat _lfeatures;
  double _loss = std::numeric_limits<double>::max(); /**< current loss. */
};

// global nn variables etc...
dMat gfeatures = dMat::Zero(784,100);
dMat glabels = dMat::Zero(10,100);
nn gmnistnn;

// objective function
FitFunc nn_of = [](const double *x, const int N)
{
  //std::copy(x,x+N,gmnistnn._allparams.begin()); // beware.
  gmnistnn._allparams.clear();
  for (int i=0;i<N;i++)
    gmnistnn._allparams.push_back(x[i]);
  gmnistnn.forward_pass(gfeatures,glabels);

  //debug
  //std::cout << "net loss= " << gmnistnn._loss << std::endl;
  //debug
  
  return gmnistnn._loss;
};

DEFINE_string(fdata,"train.csv","name of the file that contains the training data for MNIST");
DEFINE_int32(n,100,"max number of examples to train from");
DEFINE_int32(maxsolveiter,-1,"max number of optimization iterations");
DEFINE_string(fplot,"","output file for optimization log");

//TODO: train with batches.
int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  /*google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr=1;
    google::SetLogDestination(google::INFO,"");*/
  int err = load_mnist_dataset(FLAGS_fdata,FLAGS_n,gfeatures,glabels);
  if (err)
    {
      std::cout << "error loading dataset " << FLAGS_fdata << std::endl;
      exit(1);
    }
  
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  //double minloss = 1e-3; //TODO.
  //int lambda = 1e3;
  std::vector<int> lsizes = {784, 100, 10};
  gmnistnn = nn(lsizes);

  // training.
  double sigma = 2.0;
  std::vector<double> x0(gmnistnn._allparams_dim,-std::numeric_limits<double>::max());
  CMAParameters<> cmaparams(gmnistnn._allparams_dim,&x0.front(),sigma);//,lambda);
  cmaparams.set_max_iter(FLAGS_maxsolveiter);
  cmaparams._fplot = FLAGS_fplot;
  cmaparams._algo = sepCMAES;
  CMASolutions cmasols = cmaes<>(nn_of,cmaparams);
  std::cout << "status: " << cmasols._run_status << std::endl;

  // testing.
  dMat cmat = dMat::Zero(10,10);
  gmnistnn.forward_pass(gfeatures,glabels);
  //std::cout << gmnistnn._lfeatures.cols() << " / " << gmnistnn._lfeatures.rows() << std::endl;
  for (int i=0;i<gmnistnn._lfeatures.cols();i++)
    {
      dMat::Index maxv[2];
      gmnistnn._lfeatures.col(i).maxCoeff(&maxv[0]);
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
  
  //dMat epsilon = dMat::Constant(10,1,1e-20);
  double precision = diago.transpose().cwiseQuotient(col_sums).sum() / 10.0;
  double recall = diago.cwiseQuotient(row_sums).sum() / 10.0;
  /*if (std::isnan(precision))
    precision = 1.0;
  if (std::isnan(recall))
  recall = 1.0;*/
  double accuracy = diago.sum() / cmat.sum();
  double f1 = (2 * precision * recall) / (precision + recall);

  std::cout << "precision=" << precision << " / recall=" << recall << std::endl;
  std::cout << "accuracy=" << accuracy << std::endl;
  std::cout << "f1=" << f1 << std::endl;
}
