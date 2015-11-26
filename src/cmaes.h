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

#ifndef CMAES_H
#define CMAES_H

#include "esoptimizer.h"
#include "cmastrategy.h"
#include "ipopcmastrategy.h"
#include "bipopcmastrategy.h"

namespace cma = libcmaes;

namespace libcmaes
{
  template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  CMASolutions cmaes(FitFunc &func,
		     CMAParameters<TGenoPheno> &parameters,
		     ProgressFunc<CMAParameters<TGenoPheno>,CMASolutions> &pfunc=CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultPFunc,
		     GradFunc gfunc=nullptr,
		     const CMASolutions &solutions=CMASolutions(),
		     PlotFunc<CMAParameters<TGenoPheno>,CMASolutions> &pffunc=CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultFPFunc)
    {
      switch(parameters.get_algo())
	{
	case CMAES_DEFAULT:
	{
	  // if / else as object cannot be reused for now, and we're avoiding pointers.
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> cmaes_vanilla(func,parameters);
	      if (gfunc != nullptr)
		cmaes_vanilla.set_gradient_func(gfunc);
	      cmaes_vanilla.set_progress_func(pfunc);
	      cmaes_vanilla.set_plot_func(pffunc);
	      cmaes_vanilla.optimize();
	      return cmaes_vanilla.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> cmaes_vanilla(func,parameters,solutions);
	      if (gfunc != nullptr)
		cmaes_vanilla.set_gradient_func(gfunc);
	      cmaes_vanilla.set_progress_func(pfunc);
	      cmaes_vanilla.set_plot_func(pffunc);
	      cmaes_vanilla.optimize();
	      return cmaes_vanilla.get_solutions();
	    }
	}
	case IPOP_CMAES:
	{
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.set_plot_func(pffunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.set_plot_func(pffunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	}
	case BIPOP_CMAES:
	{
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.set_plot_func(pffunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.set_plot_func(pffunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	}
	case aCMAES:
	{
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> acmaes(func,parameters);
	      if (gfunc != nullptr)
		acmaes.set_gradient_func(gfunc);
	      acmaes.set_progress_func(pfunc);
	      acmaes.set_plot_func(pffunc);
	      acmaes.optimize();
	      return acmaes.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> acmaes(func,parameters,solutions);
	      if (gfunc != nullptr)
		acmaes.set_gradient_func(gfunc);
	      acmaes.set_progress_func(pfunc);
	      acmaes.set_plot_func(pffunc);
	      acmaes.optimize();
	      return acmaes.get_solutions();
	    }
	}
	case aIPOP_CMAES:
	{
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> aipop(func,parameters);
	      if (gfunc != nullptr)
		aipop.set_gradient_func(gfunc);
	      aipop.set_progress_func(pfunc);
	      aipop.set_plot_func(pffunc);
	      aipop.optimize();
	      return aipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> aipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		aipop.set_gradient_func(gfunc);
	      aipop.set_progress_func(pfunc);
	      aipop.set_plot_func(pffunc);
	      aipop.optimize();
	      return aipop.get_solutions();
	    }
	}
	case aBIPOP_CMAES:
	{
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> abipop(func,parameters);
	      if (gfunc != nullptr)
		abipop.set_gradient_func(gfunc);
	      abipop.set_progress_func(pfunc);
	      abipop.set_plot_func(pffunc);
	      abipop.optimize();
	      return abipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> abipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		abipop.set_gradient_func(gfunc);
	      abipop.set_progress_func(pfunc);
	      abipop.set_plot_func(pffunc);
	      abipop.optimize();
	      return abipop.get_solutions();
	    }
	}
	case sepCMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters);
	      if (gfunc != nullptr)
		sepcmaes.set_gradient_func(gfunc);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.set_plot_func(pffunc);
	      sepcmaes.optimize();
	      return sepcmaes.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters,solutions);
	      if (gfunc != nullptr)
		sepcmaes.set_gradient_func(gfunc);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.set_plot_func(pffunc);
	      sepcmaes.optimize();
	      return sepcmaes.get_solutions();
	    }
	}
	case sepIPOP_CMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.set_plot_func(pffunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.set_plot_func(pffunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	}
	case sepBIPOP_CMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.set_plot_func(pffunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	}
	case sepaCMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size()==0)
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters);
	      if (gfunc != nullptr)
		sepcmaes.set_gradient_func(gfunc);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.set_plot_func(pffunc);
	      sepcmaes.optimize();
	      return sepcmaes.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<CMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> sepcmaes(func,parameters,solutions);
	      if (gfunc != nullptr)
		sepcmaes.set_gradient_func(gfunc);
	      sepcmaes.set_progress_func(pfunc);
	      sepcmaes.optimize();
	      return sepcmaes.get_solutions();
	    }
	}
	case sepaIPOP_CMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size() == 0)
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.set_plot_func(pffunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> ipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		ipop.set_gradient_func(gfunc);
	      ipop.set_progress_func(pfunc);
	      ipop.optimize();
	      return ipop.get_solutions();
	    }
	}
	case sepaBIPOP_CMAES:
	{
	  if (!parameters.is_sep())
	    parameters.set_sep();
	  if (solutions.cov().size() == 0)
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.set_plot_func(pffunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	  else
	    {
	      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,TGenoPheno>,CMAParameters<TGenoPheno>,CMASolutions> bipop(func,parameters,solutions);
	      if (gfunc != nullptr)
		bipop.set_gradient_func(gfunc);
	      bipop.set_progress_func(pfunc);
	      bipop.optimize();
	      return bipop.get_solutions();
	    }
	}
	case VD_CMAES:
	{
	  if (!parameters.is_vd())
	    parameters.set_vd();
	  ESOptimizer<CMAStrategy<VDCMAUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> vdcma(func,parameters);
	  if (gfunc != nullptr)
	    vdcma.set_gradient_func(gfunc);
	  vdcma.set_progress_func(pfunc);
	  vdcma.set_plot_func(pffunc);
	  vdcma.optimize();
	  return vdcma.get_solutions();
	}
	case VD_IPOP_CMAES:
	{
	  if (!parameters.is_vd())
	    parameters.set_vd();
	  ESOptimizer<IPOPCMAStrategy<VDCMAUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> ipop(func,parameters);
	  if (gfunc != nullptr)
	    ipop.set_gradient_func(gfunc);
	  ipop.set_progress_func(pfunc);
	  ipop.set_plot_func(pffunc);
	  ipop.optimize();
	  return ipop.get_solutions();
	}
	case VD_BIPOP_CMAES:
	{
	  if (!parameters.is_vd())
	    parameters.set_vd();
	  ESOptimizer<BIPOPCMAStrategy<VDCMAUpdate,TGenoPheno>,CMAParameters<TGenoPheno>> bipop(func,parameters);
	  if (gfunc != nullptr)
	    bipop.set_gradient_func(gfunc);
	  bipop.set_progress_func(pfunc);
	  bipop.set_plot_func(pffunc);
	  bipop.optimize();
	  return bipop.get_solutions();
	}
	default:
	return CMASolutions();
	}
    }
}

#endif
