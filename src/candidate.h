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

#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "eo_matrix.h"
#include "cmaparameters.h"

namespace libcmaes
{
  /**
   * \brief candidate solution point, in function parameter space.
   */
  class Candidate
  {
  public:
    /**
     * \brief empty constructor.
     */
  Candidate():
    _fvalue(std::numeric_limits<double>::quiet_NaN()) {}
    
    /**
     * \brief constructor.
     * @param fvalue function value
     * @param x function parameter vector
     */
  Candidate(const double &fvalue,
	    const dVec &x)
    :_fvalue(fvalue),_x(x)
    {}

  ~Candidate() {}

  /**
   * \brief set candidate's function value.
   * @param fval function value
   */
  inline void set_fvalue(const double &fval) { _fvalue = fval; }

  /**
   * \brief get function value of this candidate.
   * @return function value
   */
  inline double get_fvalue() const { return _fvalue; }

  /**
   * \brief sets parameter vector of this candidate.
   * @param x parameter vector
   */
  inline void set_x(const dVec &x) { _x = x; }
  
  /**
   * \brief get parameter vector of this candidate in Eigen vector format.
   * @return parameter vector in Eigen vector format
   */
  inline dVec get_x_dvec() const { return _x; }

  /**
   * \brief get reference parameter vector of this candidate in Eigen vector format.
   * @return reference to parameter vector in Eigen vector format
   */
  inline dVec& get_x_dvec_ref() { return _x; }
  
  /**
   * \brief get parameter vector pointer of this candidate as array. 
   *        DO NOT USE from temporary candidate object.
   * @return parameter vector pointer
   */
  inline const double* get_x_ptr() const { return _x.data(); }
  
  /**
   * \brief get parameter vector copy for this candidate.
   * @return parameter vector copy
   */
  inline std::vector<double> get_x() const
  {
    std::vector<double> x;
    x.assign(_x.data(),_x.data()+_x.size());
    return x;
  }

  /**
   * \brief get x vector size
   * @return x vector size
   */
  inline unsigned int get_x_size() const { return _x.size(); }
  
  /**
   * \brief get pheno transform of parameter vector of this candidate in Eigen vector format.
   * @return pheno transform of parameter vector in Eigen vector format
   */
  template<class TGenoPheno>
    dVec get_x_pheno_dvec(const CMAParameters<TGenoPheno> &p) const
    {
      dVec gx = p.get_gp().pheno(_x);
      return gx;
    }
  
  /**
   * \brief set candidate id
   * @param id candidate id
   */
  inline void set_id(const int &id) { _id = id; }
  
  /**
   * \brief get candidate id
   * @return candidate id
   */
  inline int get_id() const { return _id; }
  
  /**
   * \brief set candidate rank
   * @param r candidate rank
   */
  inline void set_rank(const int &r) { _r = r; }

  /**
   * \brief get candidate rank
   * @return candidate rank
   */
  inline int get_rank() const { return _r; }

  protected:
   double _fvalue; /**< function value. */
   dVec _x; /**< function parameter vector. */
   int _id = -1; /**< candidate id, used for identification after ranking, when needed. */
   int _r = -1; /**< candidate rank. */
  };

  class RankedCandidate : public Candidate
  {
  public:
    RankedCandidate(const double &fvalue_mut,
		    Candidate &c,
		    const int &idx)
      :Candidate(c.get_fvalue(),dVec()),_idx(idx),_fvalue_mut(fvalue_mut)
    {}
    ~RankedCandidate() {}

    int _idx = -1;
    double _fvalue_mut;
    int _r1 = 0;
    int _r2 = 0;
    double _delta = 0.0;
  };

}

#endif
