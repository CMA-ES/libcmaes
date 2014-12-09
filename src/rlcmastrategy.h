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

#ifndef RLCMASTRATEGY_H
#define RLCMASTRATEGY_H

#include "cmastrategy.h"
#include <memory>
#include <unordered_map>
#include <iostream>

namespace libcmaes
{
  template <class T>
    class multi_vec
    {
    public:
      using param=std::vector<size_t>;

      multi_vec() {}
      explicit multi_vec(const param& dimensions)
        : dim{dimensions}, prod {1}
      {
	std::for_each(dim.begin(), dim.end(), [this] (std::size_t val)
		      {
			//std::cout << "prod=" << prod << " -- val=" << val << std::endl;
			mult.emplace_back(prod);
			prod *= val;
		      } );
        ptr.reset(new T[prod]); //TODO: deal with memory overflow
	for (int i=0;i<static_cast<int>(prod);i++)
	  ptr[i] = 1e4;
      }
      /*multi_vec(multi_vec &&mv)
      :dim(std::move(mv.dim)),mult(std::move(mv.mult)),
	prod(std::move(prod)),ptr(std::move(ptr))
	{
	}*/
      /*multi_vec(multi_vec &mv)
	:dim(mv.dim),mult(mv.mult),
	prod(mv.prod),ptr(mv.ptr)
	{
	}*/

      std::size_t capacity() const { return prod; }

      // undefined if elements in lookup != elemenets in dim
      // undefined if any element in lookup
      // is greater than or equal to corresponding dim element
      T& operator[](const param& lookup)
	{
	  return ptr[get_offset(lookup)];
	}
      const T operator[](const param& lookup) const
      {
        return ptr[get_offset(lookup)];
      }

      std::ostream& print(std::ostream &out, const param &lookup) const
	{
	  out << ptr[get_offset(lookup)];
	  return out;
	}

      std::size_t get_offset(const param& lookup) const
	{
	  std::size_t offset=0;
	  auto mit=mult.begin();
	  std::for_each(lookup.begin(), lookup.end(), [&offset, &mit] (std::size_t val)
			{
			  offset+=*mit * val;
			  ++mit;
			} );
	  return offset;
	}
      param dim;
      param mult;
      std::size_t prod;
      std::unique_ptr<T[]> ptr;
    };

  class ApproxTab
  {
  public:
    ApproxTab() {}
    ApproxTab(const int &na,
	      const int &dim,
	      const double &lbound=-5.0, const double &ubound=5.0,
	      const double &precision=1e-1);
    ~ApproxTab();
    
    double get(const dVec &s, const int &a) const;
    void set(const dVec &s, const int &a, const double &val);

    double geth(const std::string &s, const int &a) const;
    
    std::vector<size_t> vec2lookup(const dVec &s) const;
    std::string hashl(const dVec &s) const;
    
    std::ostream& print(std::ostream &out) const;
    std::ostream& print_best(std::ostream &out) const;

    std::vector<multi_vec<double>*> _tab;
    std::vector<std::unordered_map<std::string,double>> _htab;
    double _precision = 1e-1; /**< precison on tiles. */
    double _lbound = -5.0;
    double _ubound = 5.0;
    int _dim;
    int _disc; /**< discretization step. */
    int _default_val = 1e6;
    bool _use_hash = true;
  };

  class RL
  {
  public:
    RL() {}
    RL(const int &na,
       const int &dim,
       const double &lbound=-5.0, const double &ubound=5.0,
       const double &precision=1e-1);
    ~RL() 
      { 
	if (_Q)
	  delete _Q; 
      }

    int choose_action_eps_greedy(const dVec &s);
    int choose_best_action(const dVec &s);

    void max_qas(const dVec &s,
		 int &min_a, 
		 double &max_val);

    void qlearn(const dVec &s, 
		const int &a,
		const dVec &sp,
		const double &fitness);
    
    ApproxTab *_Q = nullptr; /**< Q-values. */
    double _alpha = 0.1;
    double _gamma = 0.9;
    double _epsilon = 0.1;

    std::uniform_real_distribution<double> _unif;
    std::uniform_int_distribution<> _unifi;
    std::mt19937 _gen;
  };
  
  template<class TCovarianceUpdate=CovarianceUpdate, class TGenoPheno=GenoPheno<NoBoundStrategy>>
    class RLCMAStrategy : public CMAStrategy<TCovarianceUpdate,TGenoPheno>
    {
    public:
      RLCMAStrategy(FitFunc &func,
		    CMAParameters<TGenoPheno> &parameters);
    ~RLCMAStrategy();

      void set_rl(const int &na,
		  const int &dim,
		  const double &alpha=0.1,
		  const double &gamma=0.1,
		  const double &lbound=-5.0, const double &ubound=5.0,
		  const double &precision=1e-1)
      {
	if (_rl)
	  delete _rl;
	_rl = new RL(na,dim,lbound,ubound,precision);
	_rl->_alpha = alpha;
	_rl->_gamma = gamma;
      }

      int optimize();

      void rllearn(const int &episodes); //TODO.

    void test(const int &ncontrol,
	      double &mean_fevals,
	      double &mean_reward);
      
      RL *_rl = nullptr;
    int _test_step = 100;
    };

}

#endif
