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

#include "rlcmastrategy.h"
#include "opti_err.h"
#include "llogging.h"
#include <iostream>
#include <assert.h>

namespace libcmaes
{
  /*- ApproxTab -*/
  ApproxTab::ApproxTab(const int &na,
		       const int &dim,
		       const double &lbound,
		       const double &ubound,
		       const double &precision)
    :_precision(precision),_lbound(lbound),_ubound(ubound),_dim(dim)
  {
    _disc = std::floor(ubound-lbound)/precision;
    std::cout << "disc=" << _disc + 2 << std::endl;
    for (int i=0;i<na;i++)
      {
	if (!_use_hash)
	  {
	    std::vector<size_t> dims;
	    for (int j=0;j<dim;j++)
	      dims.push_back(_disc + 2); // + 2 if for open ended interval on both sides
	    multi_vec<double> *mv = new multi_vec<double>(dims);
	    _tab.push_back(mv);
	  }
	else
	  {
	    _htab.emplace_back();
	  }
      }
  }

  ApproxTab::~ApproxTab()
  {
    if (!_use_hash)
      for (size_t i=0;i<_tab.size();i++)
	delete _tab.at(i);
  }

  double ApproxTab::get(const dVec &s, const int &a) const
  {
    std::vector<size_t> lookup = vec2lookup(s);
    if (!_use_hash)
      return (*_tab.at(a))[lookup];
    else
      {
	std::string hs = hashl(s);
	return geth(hs,a);
      }
  }
  
  void ApproxTab::set(const dVec &s, const int &a, const double &val)
  {

    if (!_use_hash)
      {
	std::vector<size_t> lookup = vec2lookup(s);
	(*_tab.at(a))[lookup] = val;
      }
    else
      {
	std::string hs = hashl(s);
	std::unordered_map<std::string,std::pair<double,int>>::iterator hit;
	if ((hit=_htab.at(a).find(hs))!=_htab.at(a).end())
	  {
	    (*hit).second.first = val;
	    (*hit).second.second++;
	  }
	else _htab.at(a).insert(std::pair<std::string,std::pair<double,int>>(hs,std::pair<double,int>(val,1)));
      }
  }

  double ApproxTab::geth(const std::string &s, const int &a) const
  {
    std::unordered_map<std::string,std::pair<double,int>>::const_iterator hit;
    if ((hit=_htab.at(a).find(s))!=_htab.at(a).end())
      {
	return (*hit).second.first;
      }
    else return -_default_val;
  }

  int ApproxTab::gethcount(const dVec &s, const int &a) const
  {
    std::string hs = hashl(s);
    std::unordered_map<std::string,std::pair<double,int>>::const_iterator hit;
    if ((hit=_htab.at(a).find(hs))!=_htab.at(a).end())
      {
	return (*hit).second.second;
      }
    else return 1;
  }
  
  std::vector<size_t> ApproxTab::vec2lookup(const dVec &s) const
  {
    std::vector<size_t> lookup;
    for (int i=0;i<s.size();i++)
      {
	double a = s[i];
	int pos = std::floor(a - _lbound) / _precision;
	if (pos < 0)
	  pos = -1;
	else if (pos > _disc)
	  pos = _disc;
	//std::cout << "pos=" << pos << " / a=" << a << " / lbound=" << _lbound << " / precision=" << _precision << std::endl;
	lookup.push_back(pos+1);
      }
    
    //debug
    /*std::cout << "s=" << s.transpose() << std::endl;
    std::cout << "lookup=";
    std::copy(lookup.begin(),lookup.end(),std::ostream_iterator<size_t>(std::cout," "));
    std::cout << std::endl;*/
    //debug

    return lookup;
  }

  std::string ApproxTab::hashl(const dVec &s) const
  {
    std::vector<size_t> lookup = vec2lookup(s);
    std::string hs;
    for (size_t c : lookup)
      {
	hs += std::to_string(c) + "|";
      }
    //std::cout << "hs=" << hs << std::endl;
    return hs;
  }

  std::ostream& ApproxTab::print(std::ostream &out) const
  {
    /*for (size_t a=0;a<_tab.size();a++)
      {
	out << "lambda=" << a + 2 << std::endl;
	dVec s = dVec::Constant(_dim,_lbound-1);
	while(s[_dim-1]<_ubound)
	  {
	    for (int i=0;i<_dim;i++)
	      {
		s[i] += 1;
		for (int j=0;j<_disc+2;j++)
		  {
		    out << get(s,a) << " ";
		  }
		out << std::endl;
	      }
	    out << std::endl;
	  }
	  }*/
    int a = 0;
    for (auto ht : _htab)
      {
	out << "lambda=" << a + 2 << std::endl;
	auto hit = ht.begin();
	while(hit!=ht.end())
	  {
	    out << (*hit).first << " / " << (*hit).second.first << std::endl;;
	    ++hit;
	  }
	++a;
      }
    return out;
  }

  std::ostream& ApproxTab::print_best(std::ostream &out) const
  {
    //TODO: get all keys.
    size_t all_keys_size = 0;
    for (auto ht: _htab)
      all_keys_size += ht.size();
    std::vector<std::string> all_keys;
    all_keys.reserve(all_keys_size);
    for (auto ht: _htab)
      {
	auto hit = ht.begin();
	while(hit!=ht.end())
	  {
	    all_keys.push_back((*hit).first);
	    ++hit;
	  }
      }
    out << "number of buckets=" << all_keys.size() << std::endl;

    // for each key, get max action
    std::unordered_map<std::string,std::vector<double>> best_actions;
    for (std::string key: all_keys)
      {
	int best_a = -1;
	int best_val = -_default_val * 1e6;
	for (size_t a=0;a<_htab.size();a++)
	  {
	    double asval = geth(key,a);
	    if (asval > best_val) // XXX: no handling of ties
	      {
		best_a = a;
		best_val = asval;
	      }
	  }
	std::vector<double> res; // 0: action, 1: qval, 2: fval
	res.push_back(best_a);
	res.push_back(best_val);
	best_actions.insert(std::pair<std::string,std::vector<double>>(key,res));
      }

    //TODO: print out max action per keys, order per fval ?
    auto bhit = best_actions.begin();
    while(bhit != best_actions.end())
      {
	out << (*bhit).first << " / " << (*bhit).second.at(0)+2 << " / " << (*bhit).second.at(1) << std::endl;
	++bhit;
      }
    
    return out;
  }

  /*- RL -*/
  RL::RL(const int &na,
	 const int &dim,
	 const double &lbound, const double &ubound,
	 const double &precision)
    :_Q(new ApproxTab(na,dim,lbound,ubound,precision))
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _unif = std::uniform_real_distribution<double>(0.0,1.0);
    _unifi = std::uniform_int_distribution<>(0,na-1);
  }

  int RL::choose_action_eps_greedy(const dVec &s)
  {
    //std::cout << "choose action Q precision=" << _Q->_precision << std::endl;

    // get current max action
    int max_a;
    double max_val;
    max_qas(s,max_a,max_val);

    // sample
    double p = _unif(_gen);
        
    // choose action.
    if (p < _epsilon) // random action
      {
	int ra = _unifi(_gen);
	return ra;
      }
    else return max_a;
  }

  int RL::choose_best_action(const dVec &s)
  {
    // get current max action
    int max_a;
    double max_val;
    max_qas(s,max_a,max_val);
    return max_a;
  }

  void RL::max_qas(const dVec &s,
		   int &max_a, 
		   double &max_val)
  {
    std::vector<int> ties;
    max_a = -1;
    max_val = -_Q->_default_val * 1e6;
    for (size_t i=0;i<_Q->_htab.size();i++)
      {
	double qval = _Q->get(s,i);
	if (qval == max_val)
	  ties.push_back(i);
	else if (qval > max_val) // XXX: no handling of ties
	  {
	    max_val = qval;
	    max_a = i;
	    ties.clear();
	    ties.push_back(max_a);
	  }
      }
    if (!ties.empty())
      {
	std::random_shuffle(ties.begin(),ties.end());
	max_a = (*ties.begin());
      }
    //assert(max_a >= 0);
  }

  void RL::qlearn(const dVec &s, 
		  const int &a,
		  const dVec &sp,
		  const double &fitness)
  {
    double qas = _Q->get(s,a);
    double max_qasp;
    int max_a;
    max_qas(sp,max_a,max_qasp);
    double alpha;
    if (_alpha < 0)
      alpha = _alpha/static_cast<double>(_Q->gethcount(sp,max_a));
    else alpha = _alpha;
    double upd_val = qas + alpha*(fitness + _gamma*max_qasp - qas);
    _Q->set(s,a,upd_val);
  }

  /*- RLCMAStrategy -*/
  template <class TCovarianceUpdate, class TGenoPheno>
  RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::RLCMAStrategy(FitFunc &func,
							     CMAParameters<TGenoPheno> &parameters)
    :CMAStrategy<TCovarianceUpdate,TGenoPheno>(func,parameters)
  {
  }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::~RLCMAStrategy()
  {
    if (_rl)
      delete _rl;
  }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  int RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::optimize()
    {
      //TODO apply optimal policy
      dVec s;
      
      CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters);
      this->_niter = 0;
      
      while(!this->stop())
	{
	  s = this->_solutions.xmean();

	  //TODO: choose best action
	  int a = _rl->choose_best_action(s);
	  
	  // translate action into lambda
	  int lambda = a+2; // XXX: basic.
	  //std::cout << "best action=" << lambda << std::endl;
	  this->_parameters._lambda = lambda;
	  this->_parameters.initialize_parameters();
	  this->_solutions._candidates.resize(lambda); // XXX: _max_hist and _kcand are lambda-dependent
	  this->_solutions._kcand = std::min(this->_parameters._lambda-1,static_cast<int>(1.0+ceil(0.1+this->_parameters._lambda/4.0)));
	
	  // run CMA step
	  dMat candidates = this->ask();
	  this->eval(candidates,this->_parameters.get_gp().pheno(candidates));
	  this->tell();
	  this->inc_iter();
	}
      
      if (this->_solutions.run_status() >= 0)
	return OPTI_SUCCESS;
      else return OPTI_ERR_TERMINATION;
    }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  void RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::rllearn(const int &episodes)
    {
      int episode = 0;

      std::cout << "Q function precision=" << _rl->_Q->_precision << std::endl;
      std::cout << "episodes\tmean_fevals\tmean_reward\n";
      
      // repeat for episodes below.
      while(episode < episodes)
	{
      	  //TODO: start from x0.
	  dVec s = this->_parameters.get_x0max(); // initial parameter values x0
	  
	  // restart.
	  //std::cout << "new episode\n";
	  CMAStrategy<TCovarianceUpdate,TGenoPheno>::_solutions = CMASolutions(CMAStrategy<TCovarianceUpdate,TGenoPheno>::_parameters);
	  this->_niter = 0;

	  // run action
	  while(!this->stop()) // run CMA steps
	    {
	      s = this->_solutions.xmean();
	      
	      // for each step, take action lambda
	      int a = _rl->choose_action_eps_greedy(s);

	      // translate action into lambda
	      int lambda = a+2; // XXX: basic.
	      this->_parameters._lambda = lambda;
	      this->_parameters.initialize_parameters();
	      if (this->_solutions._candidates.size() != lambda)
		this->_solutions._candidates.resize(lambda); // XXX: _max_hist and _kcand are lambda-dependent
	      this->_solutions._kcand = std::min(this->_parameters._lambda-1,static_cast<int>(1.0+ceil(0.1+this->_parameters._lambda/4.0)));
	      
	      //std::cout << "choosing lambda=" << lambda << std::endl;
	      
	      dMat candidates = this->ask();
	      this->eval(candidates,this->_parameters.get_gp().pheno(candidates));
	      this->tell();
	      
	      // gather fitness and update with sarsa or qlearning
	      dVec sp = this->_solutions.xmean();
	      double fitness = this->_func(sp.data(),sp.size()) + this->_solutions.fevals();
	      if (this->stop())
		fitness = 0.0;
	      _rl->qlearn(s,a,sp,-fitness); //TODO: fitness that takes fevals into account.
	    
	      //std::cout << "fitness=" << fitness << std::endl;
	      
	      this->inc_iter();
	    }

	  /*std::cout << "Q-table:\n";
	  _rl->_Q->print(std::cout);
	  std::cout << std::endl;*/

	  if (episode % _test_step == 0)
	    {
	      double mean_fevals,mean_reward;
	      test(100,mean_fevals,mean_reward);
	      //std::cout << "step=" << episode << " / fevals=" << mean_fevals << " / reward=" << -mean_reward << std::endl;
	      std::cout << episode << "\t" << mean_fevals << "\t" << -mean_reward << std::endl;
	    }
	  
	  episode++;
	}
      
      /*std::cout << "Q-table:\n";
      _rl->_Q->print_best(std::cout);
      std::cout << std::endl;*/
    }

    template <class TCovarianceUpdate, class TGenoPheno>
    void RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::test(const int &ncontrol,
							   double &mean_fevals,
							   double &mean_reward)
    {
      mean_fevals = mean_reward = 0.0;
      for (int i=0;i<ncontrol;i++)
	{
	  optimize();
	  mean_fevals += this->_solutions.fevals();
	  mean_reward += this->_solutions.best_candidate().get_fvalue() + this->_solutions.fevals();
	}
      mean_fevals /= static_cast<double>(ncontrol);
      mean_reward /= static_cast<double>(ncontrol);
    }

  template class RLCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  //template class RLCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
}
