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

bool isZero(double d) { return d == 0; }

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
  if (ne < trn)
    {
      features.resize(784,ne);
      labels.resize(10,ne);
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
int gsigmoid = 1;
bool gregularize = false;
double gtrainacc = 0.0;
double gtestacc = 0.0;

// testing
double testing(const dVec &x,
	       const bool &training=true,
	       const bool &printout=true)
{
  dMat cmat = dMat::Zero(10,10);
  //Candidate bcand = cmasols.best_candidate();
  gmnistnn._allparams.clear();
  std::copy(x.data(),x.data()+x.size(),std::back_inserter(gmnistnn._allparams));
  if (training)
    gmnistnn.forward_pass(gfeatures,glabels);
  else gmnistnn.forward_pass(gtfeatures,gtlabels);
  gmnistnn.clear_grad();
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
  
  dMat diago = cmat.diagonal();
  dMat col_sums = cmat.colwise().sum();
  dMat row_sums = cmat.rowwise().sum();
  double precision = diago.transpose().cwiseQuotient(col_sums).sum() / 10.0;
  double recall = diago.cwiseQuotient(row_sums).sum() / 10.0;
  double accuracy = diago.sum() / cmat.sum();
  double f1 = (2 * precision * recall) / (precision + recall);

  if (printout)
    {
      std::cerr << "cmat:" << std::endl << cmat << std::endl;
      if (training)
	std::cout << "training set:\n";
      else std::cout << "testing set:\n";
      std::cout << "precision=" << precision << " / recall=" << recall << std::endl;
      std::cout << "accuracy=" << accuracy << std::endl;
      std::cout << "f1=" << f1 << std::endl;
    }
  return accuracy;
}

std::random_device rd;
std::mt19937 ggen(rd());
std::uniform_int_distribution<> gunif(0,41999);
int gbatches = -1;
dMat gbatchfeatures;
dMat gbatchlabels;
std::map<int,double> gndropdims;
std::vector<double> gallparams;
double gl1reg, gl2reg;

// objective function
FitFunc nn_of = [](const double *x, const int N)
{
  nn mgn = nn(glsizes,gsigmoid,false,gregularize,gregularize,gl1reg,gl2reg);
  for (int i=0;i<N;i++)
    mgn._allparams.push_back(x[i]);
  if (gbatches <= 0)
    mgn.forward_pass(gfeatures,glabels);
  else
    mgn.forward_pass(gbatchfeatures,gbatchlabels); // gbatchfeatures set in progress function.
      
  //debug
  //std::cout << "net loss=" << gmnistnn._loss << std::endl;
  //debug
  
  return mgn._loss;
};

int gnpasses = 1;
int gi = 0;
bool gdrop = false;
FitFunc nn_dof = [](const double *x, const int N)
{
  nn mgn = nn(glsizes,gsigmoid,false,gregularize,gregularize,gl1reg,gl2reg);
  if (gdrop)
    {
      mgn._allparams = gallparams;
      int i = 0;
      auto mit = gndropdims.begin();
      while(mit!=gndropdims.end())
	{
	  mgn._allparams.at((*mit).first) = x[i];
	  ++mit;
	  ++i;
	}
    }
  else
    {
      for (int i=0;i<N;i++)
	mgn._allparams.push_back(x[i]);
    }

  dMat ggbatchfeatures;
  dMat ggbatchlabels;
  if (gbatches > 0)
    {
      int beg = gi*gbatches;
      int bsize = gbatches;
      if (gi == gnpasses-1)
	{
	  bsize = gfeatures.cols()-beg;
	  gi = 0;
	}
      else ++gi;
      ggbatchfeatures = gfeatures.block(0,beg,gfeatures.rows(),bsize);
      ggbatchlabels = glabels.block(0,beg,glabels.rows(),bsize);
    }

  if (gbatches < 0)
    mgn.forward_pass(gfeatures,glabels);
  else mgn.forward_pass(ggbatchfeatures,ggbatchlabels);

  //std::cout << "net loss=" << mgn._loss << std::endl;
  
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

double gmaxmaha = 0.0;
ProgressFunc<CMAParameters<>,CMASolutions> mpfunc = [](const CMAParameters<> &cmaparams, const CMASolutions &cmasols)
{
  if (!gdrop)
    {
      gtrainacc = testing(cmasols.best_candidate().get_x_dvec(),true,false);
      gtestacc = testing(cmasols.best_candidate().get_x_dvec(),false,false);
    }
  std::cout << "epoch=" << ceil(cmasols.niter() / gnpasses) << " / iter=" << cmasols.niter() << " / evals=" << cmaparams.lambda() * cmasols.niter() << " / f-value=" << cmasols.best_candidate().get_fvalue() << " / trainacc=" << gtrainacc << " / testacc=" << gtestacc << " / sigma=" << cmasols.sigma() /*<< " / maha=" << cmasols.maha() */ << " / iter=" << cmasols.elapsed_last_iter() << std::endl;
  
  /*if (gbatches > 0)
    {
      int beg = gi*gbatches;
      int bsize = gbatches;
      if (gi == gnpasses-1)
	{
	  bsize = gfeatures.cols()-beg;
	  gi = 0;
	}
      else ++gi;
      gbatchfeatures = gfeatures.block(0,beg,gfeatures.rows(),bsize);
      gbatchlabels = glabels.block(0,beg,glabels.rows(),bsize);
      }*/

  if (cmasols.maha() > gmaxmaha)
    gmaxmaha = cmasols.maha();
  if (cmasols.maha() < gmaxmaha / 10.0)
    {
      std::cout << "Mahalanobis drop to 1/10th of its max, stopping optimization\n";
      return 1;
    }
  
  return 0;
};
							
PlotFunc<CMAParameters<>,CMASolutions> mpffunc = [](const CMAParameters<> &cmaparams, const CMASolutions &cmasols, std::ofstream &fplotstream)
{
  std::string sep = " ";
  fplotstream << cmasols.best_candidate().get_fvalue() << sep << cmasols.fevals() << sep << cmasols.sigma() << sep << sqrt(cmasols.max_eigenv()/cmasols.min_eigenv()) << sep << cmasols.max_eigenv() << sep << gtrainacc << sep << gtestacc << std::endl;
  return 0;
};

DEFINE_string(fdata,"train.csv","name of the file that contains the training data for MNIST");
DEFINE_int32(n,100,"max number of examples to train from");
DEFINE_int32(maxsolveiter,-1,"max number of optimization iterations");
DEFINE_string(fplot,"","output file for optimization log");
DEFINE_bool(check_grad,false,"checks on gradient correctness via back propagation");
DEFINE_bool(with_gradient,false,"whether to use the gradient (backpropagation) along with black-box optimization");
DEFINE_double(lambda,-1,"number of offsprings at each generation");
DEFINE_double(sigma0,0.01,"initial value for step-size sigma (-1.0 for automated value)");
DEFINE_string(hlayers,"100","comma separated list of number of neurons per hidden layer");
DEFINE_int32(punit,1,"what neuronal processing units to use (0: sigmoid, 1: default is tanh, 2: relu)");
DEFINE_double(testp,0.0,"percentage of the training set used for testing");
DEFINE_int32(mbatch,-1,"size of minibatch");
DEFINE_int32(seed,0,"seed for es");
DEFINE_string(testf,"","test file if different than training file");
DEFINE_double(x0,-std::numeric_limits<double>::max(),"initial value for all components of the mean vector (-DBL_MAX for automated value)");
DEFINE_bool(nmbatch,false,"whether to use minibatches");
DEFINE_int32(nmbatch_budget,-1,"max budget when using minibatches");
DEFINE_double(nmbatch_ftarget,1e-3,"loss target when using minibatches");
DEFINE_bool(nmbatch_sim,false,"simplified output for minibatches in order to pipe to file");
DEFINE_double(nmbatch_acc,0.97,"accuracy target when using minibatches");
DEFINE_bool(nmbatch_rand,false,"whether to use random minibatches");
DEFINE_bool(drop,false,"whether to use dropout-like strategy for blackbox optimization");
DEFINE_int32(dropdim,100,"number of neurons being retained for optimization on each pass");
DEFINE_int32(maxdroppasses,-1,"max number of passes in drop mode");
DEFINE_bool(regularize,false,"whether to use regularization");
DEFINE_double(l1reg,0.0,"L1 regularization factor");
DEFINE_double(l2reg,1e-4,"L2 regularization weight");
DEFINE_bool(sgd,false,"run stochastic gradient descent, for comparison purpose");
DEFINE_int32(maxepochs,-1,"maximum minibatch epochs");
DEFINE_double(ftolerance,1e-12,"tolerance for comparing function values");

int main(int argc, char *argv[])
{
  ggen.seed(static_cast<uint64_t>(time(nullptr)));
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_check_grad)
    {
      FLAGS_n = 10;
      FLAGS_hlayers = "10";
    }
  std::vector<std::string> hlayers_str;
  std::vector<int> hlayers;
  tokenize(FLAGS_hlayers,hlayers_str,",");
  for (size_t i=0;i<hlayers_str.size();i++)
    hlayers.push_back(atoi(hlayers_str.at(i).c_str()));
  
  int load_size = FLAGS_n;
  if (FLAGS_nmbatch)
    load_size = 60000;
  int err = load_mnist_dataset(FLAGS_fdata,load_size,FLAGS_testp,gfeatures,glabels,gtfeatures,gtlabels);
  if (err)
    {
      std::cout << "error loading dataset " << FLAGS_fdata << std::endl;
      exit(1);
    }
  if (FLAGS_testf != "")
    {
      dMat ttfeatures, ttlabels; // dummy.
      int errt = load_mnist_dataset(FLAGS_testf,10000,false,gtfeatures,gtlabels,ttfeatures,ttlabels);
      if (errt)
	{
	  std::cout << "error loading test dataset " << FLAGS_testf << std::endl;
	  exit(1);
	}
    }

  if (FLAGS_mbatch > 0)
    {
      gbatches = FLAGS_mbatch;
      gbatchfeatures = dMat(gfeatures.rows(),gbatches);
      gbatchlabels = dMat(glabels.rows(),gbatches);
    }
  gunif = std::uniform_int_distribution<>(0,gfeatures.cols()-1);
  if (FLAGS_check_grad)
    {
      // we check on random features, but we keep the original labels.
      gfeatures.resize(784,FLAGS_n);
      gfeatures = dMat::Random(784,FLAGS_n);
    }
  gsigmoid = FLAGS_punit;
  gregularize = FLAGS_regularize;
  gl1reg = FLAGS_l1reg;
  gl2reg = FLAGS_l2reg;
  
  //debug
  /*std::cout << "gfeatures: " << gfeatures << std::endl;
    std::cout << "glabels: " << glabels << std::endl;*/
  //debug
  
  glsizes.push_back(784);
  for (size_t i=0;i<hlayers.size();i++)
    glsizes.push_back(hlayers.at(i));
  glsizes.push_back(10);
  gmnistnn = nn(glsizes,FLAGS_punit,FLAGS_check_grad || FLAGS_with_gradient || FLAGS_sgd,gregularize,gregularize,FLAGS_l1reg,FLAGS_l2reg);
  gmnistnn.to_array();
  
  if (FLAGS_check_grad)
    {
      if (gmnistnn.grad_check(gfeatures,glabels))
	std::cout << "Gradient check: OK\n";
      else std::cout << "Gradient check did fail\n";
      exit(1);
    }

  std::vector<int> dropdimbounds; // allows for sampling dimensions proportionally to each inter-layer connection space.
  std::vector<int> dropdimsize;
  if (FLAGS_drop)
    {
      int nddim = 0;
      for (size_t i=0;i<gmnistnn._lweights.size();i++)
	{
	  dropdimbounds.push_back(gmnistnn._lweights.at(i).size() + gmnistnn._lb.at(i).size());
	  /*int layer_sample_size = 0;
	  if (i != gmnistnn._lweights.size()-1)
	    layer_sample_size = ceil((dropdimbounds.back() / static_cast<double>(gmnistnn._allparams_dim)) * FLAGS_dropdim - (gmnistnn._lweights.back().size()+gmnistnn._lb.back().size())/3);
	  else
	    {
	      layer_sample_size = ceil(dropdimbounds.back()/3.0);
	      layer_sample_size += std::max(0,FLAGS_dropdim-nddim);//std::accumulate(dropdimsize.begin(),dropdimsize.end(),layer_sample_size)); 
	      std::cout << "layer_sample_size=" << layer_sample_size << std::endl;
	    }
	    dropdimsize.push_back(layer_sample_size); *///TODO: beware of floor leaving some weights unreachable.*/
	  dropdimsize.push_back(ceil((1.0/static_cast<double>(gmnistnn._lweights.size())) * FLAGS_dropdim));
	  nddim += dropdimsize.back();
	  std::cout << "bound=" << dropdimbounds.back() << " / dim=" << dropdimsize.back() << std::endl;
	}
      FLAGS_dropdim = nddim;
      std::cout << "updated dropdim to " << FLAGS_dropdim << std::endl;
    }
  
  // training.
  dMat ggfeatures;
  dMat gglabels;
  CMASolutions cmasols;
  int npasses = 1;
  if (FLAGS_nmbatch)
    npasses = ceil(gfeatures.cols()/static_cast<double>(FLAGS_n));
  if (FLAGS_mbatch > 0)
      npasses = gnpasses = ceil(gfeatures.cols()/static_cast<double>(FLAGS_mbatch));
  if (FLAGS_sgd)
    npasses = ceil(gfeatures.cols()/static_cast<double>(gbatches));
  std::vector<double> sigma0(npasses,FLAGS_sigma0);
  double fvalue = 100000.0;
  double acc = 0.0;
  int nevals = 0;
  int elapsed = 0;
  int elapsed_total = 0;
  bool init = false;
  bool run = true;
  int droppasses = 0;
  int epochs = 0;
  std::uniform_int_distribution<> dropunif(0,gmnistnn._allparams_dim-1);
  
  std::cout << "npasses=" << npasses << std::endl;
  std::cout << "dim=" << gmnistnn._allparams_dim << std::endl;

  std::chrono::time_point<std::chrono::system_clock> tstart = std::chrono::system_clock::now();
  while (run)
    {
      if (!FLAGS_drop && !FLAGS_sgd)
	{
	  for (int i=0;i<npasses;i++)
	    {
	      if (fvalue <= FLAGS_nmbatch_ftarget
		  || acc >= FLAGS_nmbatch_acc
		  || (FLAGS_nmbatch_budget != -1 && nevals >= FLAGS_nmbatch_budget))
		{
		  run = false;
		  break;
		}
	      std::vector<double> x0;
	      if (i == 0 && !init)
		{
		  ggfeatures = gfeatures;
		  gglabels = glabels;
		  if (FLAGS_x0 == -std::numeric_limits<double>::max())
		    {
		      x0 = gmnistnn._allparams;
		    }
		  else x0 = std::vector<double>(gmnistnn._allparams_dim,FLAGS_x0);
		  init = true;
		}
	      else
		{
		  Candidate bcand = cmasols.best_candidate();
		  x0.clear();
		  std::copy(bcand.get_x(),bcand.get_x()+bcand.get_x_size(),std::back_inserter(x0));
		  nn hgn = nn(glsizes,gsigmoid,false,gregularize,gregularize,FLAGS_l1reg,FLAGS_l2reg);
		  for (int i=0;i<(int)x0.size();i++)
		    hgn._allparams.push_back(x0[i]);
		  hgn.forward_pass(ggfeatures,gglabels);
		  fvalue = hgn._loss; // on ggfeatures and gglabels.
		}
	      
	      double trainacc=0.0,testacc = 0.0;
	      if (cmasols.sepcov().size())
		{
		  dVec bx = cmasols.best_candidate().get_x_dvec();
		  trainacc = testing(bx,true,false);
		  testacc = testing(bx,false,false);
		}
	      std::chrono::time_point<std::chrono::system_clock> tstop = std::chrono::system_clock::now();
	      elapsed_total = std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tstart).count();
	      if (!FLAGS_nmbatch_sim)
		std::cout << "pass #" << i << " / fvalue=" << fvalue << " / nevals=" << nevals << " / trainacc=" << trainacc << " / testacc=" << testacc << " / tim=" << elapsed/1000.0 << " / timt=" << elapsed_total/1000.0 << std::endl;
	      else std::cout << fvalue << "," << nevals << "," << trainacc << "," << testacc << "," << elapsed/1000.0 << "\t" << elapsed_total / 1000.0 << std::endl;
	      
	      int beg = i*FLAGS_n;
	      int bsize = FLAGS_n;
	      if (!FLAGS_nmbatch_rand)
		{
		  if (i == npasses-1)
		    bsize = ggfeatures.cols()-i*FLAGS_n;
		  gfeatures = ggfeatures.block(0,beg,ggfeatures.rows(),bsize);
		  glabels = gglabels.block(0,beg,gglabels.rows(),bsize);
		}
	      else
		{
		  gfeatures = dMat(ggfeatures.rows(),bsize);
		  glabels = dMat(gglabels.rows(),bsize);
		  for (int j=0;j<bsize;j++)
		    {
		      int u = gunif(ggen);
		      gfeatures.col(j) = ggfeatures.col(u);
		      glabels.col(j) = gglabels.col(u);
		    }
		}
	      
	      CMAParameters<> cmaparams(gmnistnn._allparams_dim,&x0.front()/*gmnistnn._allparams.front()*/,FLAGS_sigma0,FLAGS_mbatch > 0 ? gnpasses : FLAGS_lambda,FLAGS_seed);
	      cmaparams.set_max_iter(FLAGS_maxsolveiter);
	      cmaparams.set_fplot(FLAGS_fplot);
	      cmaparams.set_algo(sepaCMAES);
	      cmaparams.set_ftarget(1e-2);
	      //cmaparams.set_mt_feval(true);
	      /*if (FLAGS_mbatch)
		cmaparams.set_noisy();*/
	      if (FLAGS_nmbatch)
		cmaparams.set_quiet(false);//true);
	      if (FLAGS_nmbatch)
		{
		  if (!FLAGS_with_gradient)
		    cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,nullptr,cmasols,mpffunc);
		  else cmasols = cmaes<>(nn_of,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc,gnn,cmasols,mpffunc);
		}
	      else
		{
		  if (!FLAGS_with_gradient)
		    cmasols = cmaes<>(nn_dof,cmaparams,mpfunc,nullptr,cmasols,mpffunc); // was nn_of
		  else cmasols = cmaes<>(nn_dof,cmaparams,mpfunc,gnn,cmasols,mpffunc); // was nn_of
		}
	      
	      sigma0[i] = cmasols.sigma();
	      nevals += cmasols.fevals();
	      elapsed += cmasols.elapsed_time();
	      if (FLAGS_mbatch > 0)
		break;
	    }
	  if (!FLAGS_nmbatch)
	    break;
	}
      // end drop
      else if (FLAGS_drop)
	{
	  gdrop = FLAGS_drop;
	  gi = 0;
	  gmaxmaha = 0;
	  // randomly drop units
	  gmnistnn.to_array();
	  if (gallparams.empty())
	    gallparams = std::vector<double>(gmnistnn._allparams_dim,0.0);
	  std::map<int,double> olddropdims = gndropdims;
	  gndropdims.clear();
	  std::vector<double> vndropdims;
	  /*for (int d=0;d<FLAGS_dropdim;d++)
	    {
	      int u = dropunif(ggen);
	      gndropdims.insert(std::pair<int,double>(u,gallparams.at(u)));
	      vndropdims.push_back(gallparams.at(u));
	      }*/
	  
	  int zeroweights = std::count_if(gallparams.begin(),gallparams.end(),isZero);
	  double factor = 1.0;
	  if (zeroweights < gallparams.size())
	    {
	      factor = 2.0;
	      std::map<double,int,std::greater<double>> bdims;
	      for (int v=0;v<cmasols.eigenvalues().size();v++)
		bdims.insert(std::pair<double,int>(cmasols.eigenvalues()(v),v));
	      int ulim = ceil(FLAGS_dropdim / factor);
	      int p = 0;
	      std::map<double,int,std::greater<double>>::const_iterator bdmit = bdims.begin();
	      
	      while(bdmit!=bdims.end())
		{
		  int u = (*bdmit).second;
		  std::map<int,double>::const_iterator odmit = olddropdims.begin();
		  int c = 0;
		  while(odmit != olddropdims.end())
		    {
		      if (c == u)
			{
			  u = (*odmit).first;
			  break;
			}
		      ++c;
		      ++odmit;
		    }
		  gndropdims.insert(std::pair<int,double>(u,gallparams.at(u)));
		  vndropdims.push_back(gallparams.at(u));
		  ++bdmit;
		  ++p;
		  if (p >= ulim)
		    break;
		}
	      std::cout << "kept " << gndropdims.size() << " dimensions\n";
	    }
	      
	  int pos = 0;
	  //std::cout << "selected dims=";
	  for (size_t l=0;l<dropdimsize.size();l++)
	    {
	      //int beg = (l == 0) ? 0 : dropdimbounds.at(l-1);
	      dropunif = std::uniform_int_distribution<>(pos,pos+dropdimbounds.at(l)-1);
	      //std::cout << "pos=" << pos << " / bound=" << pos+dropdimbounds.at(l)-1 << std::endl;
	      pos += dropdimbounds.at(l);
	      int ulim = ceil(dropdimsize.at(l) / factor);
	      //std::cout << "ulim=" << ulim << std::endl;
	      for (int d=0;d<ulim;d++)
		{
		  int u = dropunif(ggen);
		  //std::cout << u << " ";
		  gndropdims.insert(std::pair<int,double>(u,gallparams.at(u)));
		  vndropdims.push_back(gallparams.at(u));
		}
	    }
	  //std::cout << std::endl;

	  std::cout << "zero weights=" << std::count_if(gallparams.begin(),gallparams.end(),isZero) << std::endl;
	  
	  // optimize network until convergence or maxiter etc...
	  double lambda = gnpasses;
	  CMAParameters<> cmaparams(FLAGS_dropdim,&vndropdims.front(),FLAGS_sigma0,lambda,FLAGS_seed);
	  cmaparams.set_max_iter(FLAGS_maxsolveiter);
	  cmaparams.set_fplot(FLAGS_fplot);
	  cmaparams.set_algo(aCMAES);
	  cmaparams.set_ftarget(1e-2);
	  cmaparams.set_mt_feval(true);
	  cmaparams.set_quiet(false);
	  cmaparams.set_kl(true);
	  cmaparams.set_ftolerance(FLAGS_ftolerance);
	  cmasols = cmaes<>(nn_dof,cmaparams,mpfunc);
	  nevals += cmasols.fevals();
	  
	  // update state.
	  int p = 0;
	  auto mit = gndropdims.begin();
	  while(mit!=gndropdims.end())
	    {
	      gallparams.at((*mit).first) = cmasols.best_candidate().get_x_dvec()(p);
	      ++mit;
	      ++p;
	    }
	  
	  // loop / break.
	  //double acc = 0.0;
	  std::cout << "iter=" << droppasses + 1 << " / loss=" << cmasols.best_candidate().get_fvalue() << " / fevals=" << nevals;
	  if (FLAGS_testp || FLAGS_testf != "")
	    {
	      dVec bx = Map<dVec>(&gallparams.front(),gallparams.size());
	      gtrainacc = testing(bx,true,false);
	      gtestacc = testing(bx,false,false);
	      std::cout << " / trainacc=" << gtrainacc << " / testacc=" << gtestacc;
	    }
	  std::cout << std::endl;
	  if (FLAGS_maxdroppasses > 0 && ++droppasses >= FLAGS_maxdroppasses)
	    break;
	} // end drop
      else if (FLAGS_sgd)
	{
	  for (int i=0;i<npasses;i++)
	    {
	      if (gbatches <= 0)
		{
		  std::cout << "Error: need to specify minibatch size with -mbatch\n";
		  exit(1);
		}
	      /*for (int j=0;j<gbatches;j++)
		{
		  int u = gunif(ggen);
		  gbatchfeatures.col(i*gbatches+j) = gfeatures.col(u);
		  gbatchlabels.col(j) = glabels.col(u);
		  }*/
	      int beg = i*gbatches;
	      int bsize = gbatches;
	      if (i == npasses-1)
		bsize = gfeatures.cols()-i*gbatches;
	      gbatchfeatures = gfeatures.block(0,beg,gfeatures.rows(),bsize);
	      gbatchlabels = glabels.block(0,beg,glabels.rows(),bsize);
	      gmnistnn.sgd(gbatchfeatures,gbatchlabels,0.01);
	      gmnistnn.clear_grad();
	      if (i == npasses - 1)
		{
		  dVec x = Map<dVec>(&gmnistnn._allparams.front(),gmnistnn._allparams.size());;
		  gtrainacc = testing(x,true,false);
		  gtestacc = testing(x,false,false);
		  std::cout << "epochs=" << epochs << " / iter=" << epochs * npasses + i << " / loss=" << gmnistnn._loss << " / trainacc=" << gtrainacc << " / testacc=" << gtestacc << std::endl;
		}
	    }
	  // stop after n epochs.
	  if (FLAGS_maxepochs > 0 && epochs >= FLAGS_maxepochs)
	    {
	      gallparams = gmnistnn._allparams;
	      break;
	    }
	  ++epochs;
	}
    }// end run

  if (!FLAGS_drop && !FLAGS_sgd && FLAGS_mbatch <= 0)
    {
      gfeatures = ggfeatures;
      glabels = gglabels;
    }
  
  // testing on training set.
  dVec bx;
  if (!FLAGS_drop && !FLAGS_sgd && FLAGS_mbatch <= 0)
    bx = cmasols.best_candidate().get_x_dvec();
  else bx = Map<dVec>(&gallparams.front(),gallparams.size());
  testing(bx,true);
    
  // testing on testing set, if any.
  if (FLAGS_testp || FLAGS_testf != "")
    testing(bx,false);
}
