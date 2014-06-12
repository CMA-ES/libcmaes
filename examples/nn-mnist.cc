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

double range = 1.0;
struct gen_rand {
  double factor;
public:
  gen_rand(double r=1.0) : factor(range/RAND_MAX) {}
  double operator()() {
    return rand() * factor;
  }
};

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

class nn_gradient
{
public:
  nn_gradient() {};
  ~nn_gradient() {};

  void initialize_gradients(const std::vector<int> &lsizes)
  {
    for (size_t i=0;i<lsizes.size()-1;i++)
      {
	int size_in = lsizes.at(i);
	int size_out = lsizes.at(i+1);
	dMat Gw = dMat::Zero(size_in,size_out);
	_GWs.push_back(Gw);
	dMat b = dMat::Zero(size_out,1);
	_Gbs.push_back(b);
      }
  }
  
  void accumulate(const std::vector<dMat> &GWupd,
		  const std::vector<dMat> &deltas)
  {
    int j = 0;
    for (int i=GWupd.size()-1;i>=0;i--)
      {
	_GWs.at(j) += GWupd.at(i).transpose();
	j++;
      }
    j = 0;
    for (int i=deltas.size()-1;i>=0;i--)
      {
	for (int k=0;k<deltas.at(i).cols();k++)
	  _Gbs.at(j) += deltas.at(i).col(k);
	j++;
      }
  }

  void grad_to_vec(const double &n,
		   const int &dim,
		   std::vector<double> &allgradient)
  {
    allgradient.reserve(dim);
    for (size_t i=0;i<_GWs.size();i++)
      {
	//std::cerr << "GW size=" << _GWs.at(i).size() << " / Gbs size=" << _Gbs.at(i).size() << std::endl;
	std::copy(_GWs.at(i).data(),_GWs.at(i).data()+_GWs.at(i).size(),std::back_inserter(allgradient));
	std::copy(_Gbs.at(i).data(),_Gbs.at(i).data()+_Gbs.at(i).size(),std::back_inserter(allgradient));
      }
    std::transform(allgradient.begin(),allgradient.end(),allgradient.begin(),
		   std::bind1st(std::multiplies<double>(),1.0/n));
    //std::cerr << "dim=" << dim << " / allgradient size=" << allgradient.size() << std::endl;
  }
  
  std::vector<dMat> _GWs;
  std::vector<dMat> _Gbs;
};

class nn
{
public:
  nn()
  {};
  
  nn(const std::vector<int> &lsizes,
     const bool &has_grad=false)
    :_lsizes(lsizes),_has_grad(has_grad)
  {
    for (size_t i=0;i<_lsizes.size()-1;i++)
      {
	_lweights.push_back(dMat::Random(_lsizes.at(i),_lsizes.at(i+1)));
	_lb.push_back(dVec::Zero(_lsizes.at(i+1)));
	_allparams_dim += _lweights.at(i).size() + _lsizes.at(i+1);
      }
    if (_has_grad)
      _grad.initialize_gradients(lsizes);
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

  dMat sigmoid_derivative(const dMat &X)
  {
    dMat M = sigmoid(X);
    return M.cwiseProduct((dMat::Constant(M.rows(),M.cols(),1)-M));
  }
  
  static dMat softmax(const dMat &M)
  {
    dMat maxM = M.colwise().maxCoeff();
    /*std::cerr << "M=" << M.transpose() << std::endl;
      std::cerr << "maxM=" << maxM << std::endl;*/
    dMat expM(M.rows(),M.cols());
    for (int i=0;i<M.rows();i++)
      {
	for (int j=0;j<M.cols();j++)
	  {
	    expM(i,j) = exp(M(i,j)-maxM(j));
	  }
      }
    dMat sums = expM.colwise().sum();
    //std::cerr << "sums=" << sums << std::endl;
    // div by row vector.
    for (int i=0;i<expM.cols();i++)
      {
	for (int j=0;j<expM.rows();j++)
	  {
	    expM(j,i) /= sums(0,i);
	  }
      }
    //std::cerr << "expM=" << expM << std::endl;
    return expM;
  }

  dMat get_loss(const dMat &prediction, const dMat &labels)
  {
    dMat M = labels.cwiseProduct(prediction);
    dMat Mc = M.colwise().sum();
    //std::cerr << "Mc=" << Mc << std::endl;
    dMat logM(Mc.rows(),Mc.cols());
    for (int i=0;i<Mc.rows();i++)
      {
	for (int j=0;j<Mc.cols();j++)
	  {
	    logM(i,j) = Mc(i,j) == 0 ? -20.0 : log(Mc(i,j));
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
	//std::cerr << "weights=" << _lweights.at(i).transpose() << std::endl;
	//std::cerr << "activation=" << activation.transpose() << std::endl;
	if (_has_grad)
	  _activations.push_back(activation);
	if (i != _lweights.size()-1)
	  _lfeatures = sigmoid(activation);
	else _lfeatures = softmax(activation);
	//std::cerr << "lfeatures=" << _lfeatures << std::endl;
	if (_has_grad)
	  _predicts.push_back(_lfeatures);
      }
    
    // loss.
    if (labels.size() > 0) // training mode.
      {
	if (_has_grad)
	  {
	    dMat delta = _lfeatures - labels;
	    _deltas.push_back(delta);
	  }
	/*std::cout << "features:\n";
	std::cout << lfeatures << std::endl;
	std::cout << "labels:\n";
	std::cout << labels << std::endl;
	std::cout << "delta:\n";
	std::cout << delta << std::endl;*/
	//_loss = delta.norm();
	_loss = get_loss(_lfeatures,labels).mean();
	//std::cerr << "loss=" << _loss << std::endl;
	/*if (_loss == 0.0)
	  {
	    std::cout << "prediction=" << _predicts.back() << std::endl;
	    std::cout << "M=" << labels.cwiseProduct(_lfeatures) << std::endl;
	    std::cout << "weights=" << _lweights.back() << std::endl;
	    }*/
      }
  }

  void back_propagate(const dMat &features)
  {
    std::vector<dMat> GWupd;
    int j = 0;
    for (int i=(int)_lweights.size()-1;i>=0;i--)
      {
	if (i > 0) // back propagate until hidden layer.
	  {
	    dMat delta = _lweights.at(i) * _deltas.back();
	    dMat der = sigmoid_derivative(_activations.at(i-1));
	    delta = der.cwiseProduct(delta);
	    _deltas.push_back(delta);
	  }

	dMat GW;
	if (i == 0)
	  {
	    GW = _deltas.at(j) * features.transpose();
	  }
	else
	  {
	    GW = _deltas.at(j) * _predicts.at(i-1).transpose();
	  }  
	j++;
	GWupd.push_back(GW);
      }
    _grad.accumulate(GWupd,_deltas);
  }
  
  void to_array()
  {
    _allparams.clear();
    _allparams.reserve(_allparams_dim);
    auto vit = _allparams.begin();
    for (size_t i=0;i<_lweights.size();i++)
      {
	std::copy(_lweights.at(i).data(),_lweights.at(i).data()+_lweights.at(i).size(),vit); //TODO: use std::back_inserter instead.
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

  dVec grad_to_vec(const double &n)
  {
    _allgradient.clear();
    _grad.grad_to_vec(n,_allparams.size(),_allgradient);
    dVec grad(_allparams.size());
    std::copy(_allgradient.begin(),_allgradient.end(),grad.data());
    return grad;
  }
  
  void clear_grad()
  {
    _grad = nn_gradient();
    _grad.initialize_gradients(_lsizes);
    _activations.clear();
    _deltas.clear();
    _predicts.clear();
  }

  bool grad_check(const dMat &gfeatures, const dMat &glabels)
  {
    double epsilon = 1e-4;
    int attempts = 1;
    while (attempts > 0)
      {
	_allparams.reserve(_allparams_dim);
	std::generate_n(std::back_inserter(_allparams), _allparams_dim, gen_rand());
	forward_pass(gfeatures,glabels);
	double cost = _loss;
	back_propagate(gfeatures);
	grad_to_vec(gfeatures.cols()); // in allgradient

	//std::cout << "computing numerical gradient\n";
	std::vector<double> allparams = _allparams;
	_has_grad = false;
	std::vector<double> numerical_gradient;
	for (int i=0;i<(int)_allparams_dim;i++)
	  {
	    /*if (i < 784 && gfeatures(i,0) == 0)
	      {
		numerical_gradient.push_back(0.0);
		continue;
		}*/
	    std::vector<double> e(_allparams_dim,0.0);
	    e.at(i) = 2.0*epsilon;
	    std::vector<double> y1;
	    std::transform(allparams.begin(),allparams.end(),e.begin(),std::back_inserter(y1),std::plus<double>());
	    _allparams = y1;
	    //std::cerr << "old weight=" << allparams.at(i) << " / new weight=" << y1.at(i) << std::endl;
	    forward_pass(gfeatures,glabels);
	    double cost_epsilon = _loss;
	    double res = (cost_epsilon - cost) / (2.0*epsilon);
	    //std::cerr << "#i=" << i << " / cost=" << cost << " / cost_epsilon=" << cost_epsilon << " / res=" << res << std::endl;
	    numerical_gradient.push_back(res);
	  }

	double diff=0,total_diff=0,max_diff=-1.0;
	for (int i=0;i<(int)_allparams_dim;i++)
	  {
	    diff = fabs(numerical_gradient.at(i)-_allgradient.at(i));
	    std::cout << "i: " << i << " -- diff: " << diff
		      << " -- numerical gradient: " << numerical_gradient.at(i)
		      << " -- returned gradient: " << _allgradient.at(i) << std::endl;
	    max_diff = std::max(max_diff,diff);
	    total_diff += diff;
	  }
	std::cout << "attempt #" << attempts << " -- total diff: " << total_diff
		  << " -- max diff: " << max_diff << std::endl;
	
	if (max_diff > epsilon)
	  {
	    std::cout << "gradient cal failed: max diff too high\n";
	    return false;
	  }
	
	--attempts;
      }
    return true;
  }
  
  std::vector<int> _lsizes; /**< layer sizes. */
  std::vector<dMat> _lweights; /**< weight matrice, per layer. */
  std::vector<dVec> _lb;
  unsigned int _allparams_dim = 0;
  std::vector<double> _allparams; /**< all parameters, flat representation. */
  dMat _lfeatures;
  double _loss = std::numeric_limits<double>::max(); /**< current loss. */
  nn_gradient _grad;
  bool _has_grad = false;
  std::vector<double> _allgradient; /**< flat representation. */
  std::vector<dMat> _activations; // only for bp.
  std::vector<dMat> _deltas;
  std::vector<dMat> _predicts;
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
  if (gmnistnn._has_grad)
    gmnistnn.clear_grad();
  for (int i=0;i<N;i++)
    gmnistnn._allparams.push_back(x[i]);
  gmnistnn.forward_pass(gfeatures,glabels);

  //debug
  //std::cout << "net loss= " << gmnistnn._loss << std::endl;
  //debug
  
  return gmnistnn._loss;
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

//TODO: train with batches.
int main(int argc, char *argv[])
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  int hlayer = 100;
  if (FLAGS_check_grad)
    {
      FLAGS_n = 10;
      hlayer = 10;
    }

  int err = load_mnist_dataset(FLAGS_fdata,FLAGS_n,gfeatures,glabels);
  if (err)
    {
      std::cout << "error loading dataset " << FLAGS_fdata << std::endl;
      exit(1);
    }
  if (FLAGS_check_grad)
    {
      // we check on random features, but we keep the original labels.
      gfeatures.resize(784,FLAGS_n);
      gfeatures = dMat::Random(784,FLAGS_n);
    }
    
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  //double minloss = 1e-3; //TODO.
  //int lambda = 1e3;
  std::vector<int> lsizes = {784, hlayer, 10};
  gmnistnn = nn(lsizes,FLAGS_check_grad || FLAGS_with_gradient);

  if (FLAGS_check_grad)
    {
      if (gmnistnn.grad_check(gfeatures,glabels))
	std::cout << "Gradient check: OK\n";
      else std::cout << "Gradient check did fail\n";
      exit(1);
    }
  
  // training.
  double sigma = 2.0;
  std::vector<double> x0(gmnistnn._allparams_dim,-std::numeric_limits<double>::max());
  CMAParameters<> cmaparams(gmnistnn._allparams_dim,&x0.front(),sigma,FLAGS_lambda);
  cmaparams.set_max_iter(FLAGS_maxsolveiter);
  cmaparams._fplot = FLAGS_fplot;
  cmaparams._algo = sepaCMAES;
  cmaparams.set_ftarget(1e-8);
  CMASolutions cmasols;
  if (!FLAGS_with_gradient)
    cmasols = cmaes<>(nn_of,cmaparams);
  else cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,gnn);
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
