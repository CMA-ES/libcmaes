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
using namespace libcmaes;

#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/list.hpp>
using namespace boost::python;
#include "py_boost_function.hpp"

/*- required wrappers -*/
double fitfunc_f(const boost::python::list &x, const int &n)
{
  std::cout << "uninstanciated fitfunc_f\n";
  return 0.0;
}
boost::function<double(const boost::python::list&,const int&)> fitfunc_bf(fitfunc_f);

/* wrapper to cmaes high level function. */
template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  CMASolutions pcmaes(boost::function<double(const boost::python::list&,const int&)>& fitfunc_bf,
  CMAParameters<TGenoPheno> &parameters)
  {
    FitFunc fpython = [fitfunc_bf](const double *x, const int N)
      {
	boost::python::list plx;
	for (int i=0;i<N;i++) // XXX: how to avoid the loop ?
	  plx.append(x[i]);
	return fitfunc_bf(plx,N);
      };
    return cmaes(fpython,parameters);
  }

/* wrapper to CMAParameters constructor with vector. */
template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  CMAParameters<TGenoPheno> make_parameters(const boost::python::list &x0,
					    const double &sigma,
					    const int &lambda=-1,
					    const uint64_t &seed=0)/*,
					    const TGenoPheno &gp=GenoPheno<NoBoundStrategy>())*/
{
  std::vector<double> vx0;
  for (int i=0;i<len(x0);i++)
    vx0.push_back(boost::python::extract<double>(x0[i]));
  return CMAParameters<TGenoPheno>(vx0,sigma,lambda,seed);//,gp);
}

BOOST_PYTHON_MODULE(lcmaes)
{
  /*- parameters object and maker -*/
  def("make_parameters",make_parameters<GenoPheno<NoBoundStrategy>>,args("x0","sigma","lambda","seed"));
  class_<CMAParameters<GenoPheno<NoBoundStrategy>>>("CMAParametersNB")
    .def(init<const int&,const double*,const double&,const int&,const uint64_t&,const GenoPheno<NoBoundStrategy>&>())
    .def("initialize_parameters", &CMAParameters<GenoPheno<NoBoundStrategy>>::initialize_parameters)
    .def("set_noisy", &CMAParameters<GenoPheno<NoBoundStrategy>>::set_noisy)
    .def("set_sep",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_sep)
    .def("set_automaxiter",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_automaxiter)
    .def("set_fixed_p",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_fixed_p)
    .def("unset_fixed_p",&CMAParameters<GenoPheno<NoBoundStrategy>>::unset_fixed_p)
    .def("set_restarts",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_restarts)
    .def("set_max_iter",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_max_iter)
    .def("get_max_iter",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_max_iter)
    .def("set_max_fevals",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_max_fevals)
    .def("set_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_ftarget)
    .def("reset_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy>>::reset_ftarget)
    .def("get_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_ftarget)
    .def("set_seed",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_seed)
    .def("get_seed",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_seed)
    .def("set_ftolerance",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_ftolerance)
    .def("get_ftolerance",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_ftolerance)
    .def("set_xtolerance",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_xtolerance)
    .def("get_xtolerance",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_xtolerance)
    .def("lambda",&CMAParameters<GenoPheno<NoBoundStrategy>>::lambda)
    .def("dim",&CMAParameters<GenoPheno<NoBoundStrategy>>::dim)
    .def("set_quiet",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_quiet)
    .def("quiet",&CMAParameters<GenoPheno<NoBoundStrategy>>::quiet)
    .def("set_algo",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_algo)
    .def("get_algo",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_algo)
    //TODO: geno with set_gp / get_gp
    .def("set_fplot",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_fplot)
    .def("get_fplot",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_fplot)
    .def("set_gradient",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_gradient)
    .def("get_gradient",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_gradient)
    .def("set_edm",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_edm)
    .def("get_edm",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_edm)
    .def("set_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_mt_feval)
    .def("get_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_mt_feval)
    ;

  /*- FitFunc -*/  
  def_function<double(const boost::python::list&,const int&)>("fitfunc_pbf","fitfunc for python");
  scope().attr("fitfunc_bf") = fitfunc_bf;
  
  /*- solutions object -*/
  class_<CMASolutions>("CMASolutions")
    .def(init<Parameters<GenoPheno<NoBoundStrategy>>&>())
    .def("sort_candidates",&CMASolutions::sort_candidates)
    .def("update_best_candidate",&CMASolutions::update_best_candidates)
    //TODO: update_eigenv ?
    //TODO: best_candidate, need to return a python candidate object...
    .def("size",&CMASolutions::size)
    .def("edm",&CMASolutions::edm)
    //TODO: cov, see http://eigen.tuxfamily.org/index.php?title=PythonInteropBoost for wrapping Eigen3 structures.
    //TODO: sepcov
    .def("sigma",&CMASolutions::sigma)
    //TODO: xmean
    .def("run_status",&CMASolutions::run_status)
    .def("elapsed_time",&CMASolutions::elapsed_time)
    .def("niter",&CMASolutions::niter)
    ;

  /*- TODO: candidate object -*/
  class_<Candidate>("Candidate")
    //TODO: constructor with dVec ?
    .def("get_fvalue",&Candidate::get_fvalue)
    //TODO: various get_x !
    ;
  
  /* esoptimizer object -*/
  //class_<ESOptimizer<CMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>,CMAParameters<GenoPheno<NoBoundStrategy>>>>("ESOptimizer",init<FitFunc&,CMAParameters<GenoPheno<NoBoundStrategy>>&>())
	 //.def("optimize",&ESOptimizer::optimize)
  //;

  /*- cmaes header -*/
  def("pcmaes",pcmaes<GenoPheno<NoBoundStrategy>>,args("fitfunc","parameters"));
  
} // end boost
