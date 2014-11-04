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
#include "surrogatestrategy.h"
using namespace libcmaes;

#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/list.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
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
					    const uint64_t &seed=0,
					    const TGenoPheno &gp=GenoPheno<NoBoundStrategy>())
{
  std::vector<double> vx0;
  for (int i=0;i<len(x0);i++)
    vx0.push_back(boost::python::extract<double>(x0[i]));
  return CMAParameters<TGenoPheno>(vx0,sigma,lambda,seed,gp);
}

CMAParameters<GenoPheno<NoBoundStrategy>> make_simple_parameters(const boost::python::list &x0,
								 const double &sigma,
								 const int &lambda=-1,
								 const uint64_t &seed=0)
{
  std::vector<double> vx0;
  for (int i=0;i<len(x0);i++)
    vx0.push_back(boost::python::extract<double>(x0[i]));
  return CMAParameters<>(vx0,sigma,lambda,seed);
}

boost::python::list get_solution_xmean(const CMASolutions &s)
{
  boost::python::list xmean;
  for (int i=0;i<s.xmean().size();i++)
    xmean.append(s.xmean()(i));
  return xmean;
}

boost::python::list get_candidate_x(const Candidate &c)
{
  boost::python::list x;
  for (int i=0;i<(int)c.get_x_size();i++)
    x.append(c.get_x()[i]);
  return x;
}

PyObject* get_solution_cov_py(const CMASolutions &s)
{
  npy_intp shape[2] = {s.dim(),s.dim()};
  PyObject* pyArray = PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, (double*)s.cov_data());
  return pyArray;
}

boost::python::object get_solution_cov(const CMASolutions &s)
{
  PyObject *pobj = get_solution_cov_py(s);
  boost::python::object boostobj(boost::python::handle<>((PyObject*)pobj));
  return boostobj;
}

PyObject* get_solution_sepcov_py(const CMASolutions &s)
{
  npy_intp shape[1] = {s.dim()};
  PyObject* pyArray = PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, (double*)s.sepcov_data());
  return pyArray;
}

boost::python::object get_solution_sepcov(const CMASolutions &s)
{
  PyObject *pobj = get_solution_cov_py(s);
  boost::python::object boostobj(boost::python::handle<>((PyObject*)pobj));
  return boostobj;
}

template <class TBoundStrategy=NoBoundStrategy,class TScalingStrategy=NoScalingStrategy>
GenoPheno<TBoundStrategy,TScalingStrategy> make_genopheno(const boost::python::list &lbounds,
							  const boost::python::list &ubounds,
							  const int &dim)
{
  assert(len(lbounds)==len(ubounds));
  std::vector<double> vlbounds, vubounds;
  for (int i=0;i<len(lbounds);i++)
    {
      vlbounds.push_back(boost::python::extract<double>(lbounds[i]));
      vubounds.push_back(boost::python::extract<double>(ubounds[i]));
    }
  return GenoPheno<TBoundStrategy,TScalingStrategy>(&vlbounds.front(),&vubounds.front(),dim);
}

/* wrapper to surrogate cmaes high level function. */
int csurrfunc_f(const boost::python::list &candidates, boost::python::object &pyArray)
{
  std::cout << "uninstanciated csurrfunc_f\n";
  return 0.0;
}
boost::function<int(const boost::python::list&, boost::python::object &pyArray)> csurrfunc_bf(csurrfunc_f);

int surrfunc_f(const boost::python::list &candidates, boost::python::object &pyArray)
{
  std::cout << "uninstanciated surrfunc_f\n";
  return 0.0;
}
boost::function<int(boost::python::list&, boost::python::object &pyArray)> surrfunc_bf(surrfunc_f);

template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  CMASolutions surrpcmaes(boost::function<double(const boost::python::list&,const int&)>& fitfunc_bf,
			  boost::function<int(const boost::python::list&,boost::python::object&)> &csurrfunc_bf,
			  boost::function<int(boost::python::list&,boost::python::object&)> &surrfunc_bf,
			  CMAParameters<TGenoPheno> &parameters)
  {
    FitFunc fpython = [fitfunc_bf](const double *x, const int N)
      {
	boost::python::list plx;
	for (int i=0;i<N;i++)
	  plx.append(x[i]);
	return fitfunc_bf(plx,N);
      };
    ESOptimizer<ACMSurrogateStrategy<CMAStrategy,CovarianceUpdate>,CMAParameters<>> optim(fpython,parameters);
    CSurrFunc csurrfpython = [csurrfunc_bf](const std::vector<Candidate> &c, const dMat &m)
      {
	boost::python::list plc;
	for (size_t i=0;i<c.size();i++)
	  plc.append(c[i]);
	npy_intp shape[2] = {m.rows(),m.cols()};
	PyObject *pyArray = PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, (double*)m.data());
	boost::python::object boostobj(boost::python::handle<>((PyObject*)pyArray));
	return csurrfunc_bf(plc,boost::ref(boostobj));
      };
    optim.set_ftrain(csurrfpython);
    SurrFunc surrfpython = [surrfunc_bf](const std::vector<Candidate> &c, const dMat &m)
      {
	boost::python::list plc;
	for (size_t i=0;i<c.size();i++)
	  plc.append(c[i]);
	npy_intp shape[2] = {m.rows(),m.cols()};
	PyObject *pyArray = PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, (double*)m.data());
	boost::python::object boostobj(boost::python::handle<>((PyObject*)pyArray));
	return surrfunc_bf(boost::ref(plc),boost::ref(boostobj));
      };
      optim.set_fpredict(surrfpython);
    optim.optimize();
    return optim.get_solutions();
  }

BOOST_PYTHON_MODULE(lcmaes)
{
  import_array(); // numpy.
  
  /*- parameters object and maker -*/
  class_<CMAParameters<GenoPheno<NoBoundStrategy>>>("CMAParametersNB")
    .def("initialize_parameters", &CMAParameters<GenoPheno<NoBoundStrategy>>::initialize_parameters)
    .def("set_noisy", &CMAParameters<GenoPheno<NoBoundStrategy>>::set_noisy)
    .def("set_sep",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_sep)
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
    .def("set_str_algo",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_str_algo)
    .def("get_algo",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_algo)
    .def("set_fplot",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_fplot)
    .def("get_fplot",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_fplot)
    .def("set_gradient",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_gradient)
    .def("get_gradient",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_gradient)
    .def("set_edm",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_edm)
    .def("get_edm",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_edm)
    .def("set_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy>>::set_mt_feval)
    .def("get_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy>>::get_mt_feval)
    ;
  def("make_parameters",make_parameters<GenoPheno<NoBoundStrategy>>,args("x0","sigma","lambda","seed","gp"));
  def("make_simple_parameters",make_simple_parameters,args("x0","sigma","lambda","seed"));
  class_<CMAParameters<GenoPheno<pwqBoundStrategy>>>("CMAParametersNB")
    .def("initialize_parameters", &CMAParameters<GenoPheno<pwqBoundStrategy>>::initialize_parameters)
    .def("set_noisy", &CMAParameters<GenoPheno<pwqBoundStrategy>>::set_noisy)
    .def("set_sep",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_sep)
    .def("set_fixed_p",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_fixed_p)
    .def("unset_fixed_p",&CMAParameters<GenoPheno<pwqBoundStrategy>>::unset_fixed_p)
    .def("set_restarts",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_restarts)
    .def("set_max_iter",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_max_iter)
    .def("get_max_iter",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_max_iter)
    .def("set_max_fevals",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_max_fevals)
    .def("set_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_ftarget)
    .def("reset_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy>>::reset_ftarget)
    .def("get_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_ftarget)
    .def("set_seed",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_seed)
    .def("get_seed",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_seed)
    .def("set_ftolerance",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_ftolerance)
    .def("get_ftolerance",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_ftolerance)
    .def("set_xtolerance",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_xtolerance)
    .def("get_xtolerance",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_xtolerance)
    .def("lambda",&CMAParameters<GenoPheno<pwqBoundStrategy>>::lambda)
    .def("dim",&CMAParameters<GenoPheno<pwqBoundStrategy>>::dim)
    .def("set_quiet",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_quiet)
    .def("quiet",&CMAParameters<GenoPheno<pwqBoundStrategy>>::quiet)
    .def("set_str_algo",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_str_algo)
    .def("get_algo",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_algo)
    .def("set_fplot",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_fplot)
    .def("get_fplot",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_fplot)
    .def("set_gradient",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_gradient)
    .def("get_gradient",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_gradient)
    .def("set_edm",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_edm)
    .def("get_edm",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_edm)
    .def("set_mt_feval",&CMAParameters<GenoPheno<pwqBoundStrategy>>::set_mt_feval)
    .def("get_mt_feval",&CMAParameters<GenoPheno<pwqBoundStrategy>>::get_mt_feval)
    ;
  def("make_parameters_pwqb",make_parameters<GenoPheno<pwqBoundStrategy>>,args("x0","sigma","lambda","gp"));
  class_<CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>>("CMAParametersNB")
    .def("initialize_parameters", &CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::initialize_parameters)
    .def("set_noisy", &CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_noisy)
    .def("set_sep",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_sep)
    .def("set_fixed_p",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_fixed_p)
    .def("unset_fixed_p",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::unset_fixed_p)
    .def("set_restarts",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_restarts)
    .def("set_max_iter",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_max_iter)
    .def("get_max_iter",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_max_iter)
    .def("set_max_fevals",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_max_fevals)
    .def("set_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_ftarget)
    .def("reset_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::reset_ftarget)
    .def("get_ftarget",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_ftarget)
    .def("set_seed",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_seed)
    .def("get_seed",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_seed)
    .def("set_ftolerance",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_ftolerance)
    .def("get_ftolerance",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_ftolerance)
    .def("set_xtolerance",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_xtolerance)
    .def("get_xtolerance",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_xtolerance)
    .def("lambda",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::lambda)
    .def("dim",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::dim)
    .def("set_quiet",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_quiet)
    .def("quiet",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::quiet)
    .def("set_str_algo",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_str_algo)
    .def("get_algo",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_algo)
    .def("set_fplot",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_fplot)
    .def("get_fplot",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_fplot)
    .def("set_gradient",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_gradient)
    .def("get_gradient",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_gradient)
    .def("set_edm",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_edm)
    .def("get_edm",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_edm)
    .def("set_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::set_mt_feval)
    .def("get_mt_feval",&CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>::get_mt_feval)
    ;
  def("make_parameters_ls",make_parameters<GenoPheno<NoBoundStrategy,linScalingStrategy>>,args("x0","sigma","lambda","gp"));
  class_<CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>>("CMAParametersNB")
    .def("initialize_parameters", &CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::initialize_parameters)
    .def("set_noisy", &CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_noisy)
    .def("set_sep",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_sep)
    .def("set_fixed_p",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_fixed_p)
    .def("unset_fixed_p",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::unset_fixed_p)
    .def("set_restarts",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_restarts)
    .def("set_max_iter",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_max_iter)
    .def("get_max_iter",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_max_iter)
    .def("set_max_fevals",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_max_fevals)
    .def("set_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_ftarget)
    .def("reset_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::reset_ftarget)
    .def("get_ftarget",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_ftarget)
    .def("set_seed",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_seed)
    .def("get_seed",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_seed)
    .def("set_ftolerance",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_ftolerance)
    .def("get_ftolerance",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_ftolerance)
    .def("set_xtolerance",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_xtolerance)
    .def("get_xtolerance",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_xtolerance)
    .def("lambda",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::lambda)
    .def("dim",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::dim)
    .def("set_quiet",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_quiet)
    .def("quiet",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::quiet)
    .def("set_str_algo",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_str_algo)
    .def("get_algo",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_algo)
    .def("set_fplot",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_fplot)
    .def("get_fplot",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_fplot)
    .def("set_gradient",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_gradient)
    .def("get_gradient",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_gradient)
    .def("set_edm",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_edm)
    .def("get_edm",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_edm)
    .def("set_mt_feval",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::set_mt_feval)
    .def("get_mt_feval",&CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>::get_mt_feval)
    ;
    def("make_parameters_pwqb_ls",make_parameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>>,args("x0","sigma","lambda","gp"));
    
  /*- FitFunc -*/  
  def_function<double(const boost::python::list&,const int&)>("fitfunc_pbf","fitfunc for python");
  scope().attr("fitfunc_bf") = fitfunc_bf;
  
  /*- solutions object -*/
  class_<CMASolutions>("CMASolutions")
    .def(init<Parameters<GenoPheno<NoBoundStrategy>>&>())
    .def("sort_candidates",&CMASolutions::sort_candidates)
    .def("update_best_candidate",&CMASolutions::update_best_candidates)
    .def("best_candidate",&CMASolutions::best_candidate)
    .def("size",&CMASolutions::size)
    .def("edm",&CMASolutions::edm)
    .def("sigma",&CMASolutions::sigma)
    .def("fevals",&CMASolutions::fevals)
    .def("eigenvalues",&CMASolutions::eigenvalues)
    .def("min_eigenv",&CMASolutions::min_eigenv)
    .def("max_eigenv",&CMASolutions::max_eigenv)
    .def("run_status",&CMASolutions::run_status)
    .def("elapsed_time",&CMASolutions::elapsed_time)
    .def("elapsed_last_iter",&CMASolutions::elapsed_last_iter)
    .def("niter",&CMASolutions::niter)
    ;
  def("get_solution_xmean",get_solution_xmean,args("sol"));
  def("get_solution_cov",get_solution_cov,args("sol"));
  def("get_solution_sepcov",get_solution_sepcov,args("sol"));
  
  /*- solution candidate object -*/
  class_<Candidate>("Candidate")
    .def("get_fvalue",&Candidate::get_fvalue)
    ;
  def("get_candidate_x",get_candidate_x,args("cand"));

  /*- genopheno object -*/
  class_<GenoPheno<NoBoundStrategy>>("GenoPhenoNB")
    ;
  def("make_genopheno",make_genopheno<NoBoundStrategy,NoScalingStrategy>,args("lbounds","ubounds","dim"));
  class_<GenoPheno<pwqBoundStrategy>>("GenoPhenoPWQB")
    ;
  def("make_genopheno_pwqb",make_genopheno<pwqBoundStrategy,NoScalingStrategy>,args("lbounds","ubounds","dim"));
  class_<GenoPheno<NoBoundStrategy,linScalingStrategy>>("GenoPhenoLS")
    ;
  def("make_genopheno_ls",make_genopheno<NoBoundStrategy,linScalingStrategy>,args("lbounds","ubounds","dim"));
  class_<GenoPheno<pwqBoundStrategy,linScalingStrategy>>("GenoPhenoPWQBLS")
    ;
  def("make_genopheno_pwqb_ls",make_genopheno<pwqBoundStrategy,linScalingStrategy>,args("lbounds","ubounds","dim"));
  //TODO: geno/pheno custom transforms.

  
  /* esoptimizer object -*/
  //class_<ESOptimizer<CMAStrategy<CovarianceUpdate,GenoPheno<NoBoundStrategy>>,CMAParameters<GenoPheno<NoBoundStrategy>>>>("ESOptimizer",init<FitFunc&,CMAParameters<GenoPheno<NoBoundStrategy>>&>())
	 //.def("optimize",&ESOptimizer::optimize)
  //;

  /*- cmaes header -*/
  def("pcmaes",pcmaes<GenoPheno<NoBoundStrategy>>,args("fitfunc","parameters"));
  def("pcmaes_pwqb",pcmaes<GenoPheno<pwqBoundStrategy>>,args("fitfunc","parameters"));
  def("pcmaes_ls",pcmaes<GenoPheno<NoBoundStrategy,linScalingStrategy>>,args("fitfunc","parameters"));
  def("pcmaes_pwqb_ls",pcmaes<GenoPheno<pwqBoundStrategy,linScalingStrategy>>,args("fitfunc","parameters"));


  /*- surrogates -*/
  def_function<int(const boost::python::list&,boost::python::object&)>("csurrfunc_pbf","training function for python");
  scope().attr("csurrfunc_bf") = csurrfunc_bf;
  def_function<int(boost::python::list&,boost::python::object&)>("surrfunc_pbf","prediction function for python");
  scope().attr("surrfunc_bf") = surrfunc_bf;

  def("surrpcmaes",surrpcmaes<GenoPheno<NoBoundStrategy>>,args("fitfunc","trainfunc","predictfunc","parameters"));
  
} // end boost
