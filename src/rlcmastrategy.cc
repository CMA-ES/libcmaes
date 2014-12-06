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
#include "llogging.h"
#include <iostream>

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
    for (int i=0;i<na;i++)
      {
	std::vector<size_t> dims;
	for (int j=0;j<dim;j++)
	  dims.push_back(_disc + 2); // + 2 if for open ended interval on both sides
	multi_vec<double> *mv = new multi_vec<double>(dims);
	_tab.push_back(mv);
      }
  }

  ApproxTab::~ApproxTab()
  {
    for (size_t i=0;i<_tab.size();i++)
      delete _tab.at(i);
  }

  double ApproxTab::get(const dVec &s, const int &a) const
  {
    std::vector<size_t> lookup = vec2lookup(s);
    return (*_tab.at(a))[lookup];
  }
  
  void ApproxTab::set(const dVec &s, const int &a, const double &val)
  {
    std::vector<size_t> lookup = vec2lookup(s);
    (*_tab.at(a))[lookup] = val;
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

  std::ostream& ApproxTab::print(std::ostream &out) const
  {
    for (size_t a=0;a<_tab.size();a++)
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

    // get current min action
    int min_a;
    double min_val;
    min_qas(s,min_a,min_val);

    // sample
    double p = _unif(_gen);
        
    // choose action.
    if (p < _epsilon) // random action
      {
	int ra = _unifi(_gen);
	return ra;
      }
    else return min_a;
  }

  void RL::min_qas(const dVec &s,
		   int &min_a, 
		   double &min_val)
  {
    min_a = -1;
    min_val = 1e10;
    for (size_t i=0;i<_Q->_tab.size();i++)
      {
	double qval = _Q->get(s,i);
	if (qval < min_val)
	  {
	    min_val = qval;
	    min_a = i;
	  }
      }
  }

  void RL::qlearn(const dVec &s, 
		  const int &a,
		  const dVec &sp,
		  const double &fitness)
  {
    double qas = _Q->get(s,a);
    double min_qasp;
    int min_a;
    min_qas(sp,min_a,min_qasp);
    double upd_val = qas + _alpha*(fitness + _gamma*min_qasp - qas);
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
      //TODO: learn optimal policy
      rllearn();

      //TODO apply optimal policy

      return 0;
    }
  
  template <class TCovarianceUpdate, class TGenoPheno>
  void RLCMAStrategy<TCovarianceUpdate,TGenoPheno>::rllearn()
    {
      int episode = 0;

      std::cout << "Q function precision=" << _rl->_Q->_precision << std::endl;
      
      // repeat for episodes below.
      while(episode < _max_episodes)
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
	      this->_solutions._candidates.resize(lambda); // XXX: _max_hist and _kcand are lambda-dependent
	      this->_solutions._kcand = std::min(this->_parameters._lambda-1,static_cast<int>(1.0+ceil(0.1+this->_parameters._lambda/4.0)));
	      
	      //std::cout << "choosing lambda=" << lambda << std::endl;
	      
	      dMat candidates = this->ask();
	      this->eval(candidates,this->_parameters.get_gp().pheno(candidates));
	      this->tell();
	      
	      // gather fitness and update with sarsa or qlearning
	      dVec sp = this->_solutions.xmean();
	      double fitness = this->_func(sp.data(),sp.size()) + this->_solutions.fevals();
	      _rl->qlearn(s,a,sp,fitness); //TODO: fitness that takes fevals into account.
	    
	      //std::cout << "fitness=" << fitness << std::endl;
	      
	      this->inc_iter();
	    }

	  /*std::cout << "Q-table:\n";
	  _rl->_Q->print(std::cout);
	  std::cout << std::endl;*/
	  
	  episode++;
	}

      std::cout << "Q-table:\n";
      _rl->_Q->print(std::cout);
      std::cout << std::endl;
    }

  template class RLCMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>;
  template class RLCMAStrategy<CovarianceUpdate,GenoPheno<pwqBoundStrategy>>;
}
