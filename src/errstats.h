/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 INRIA
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

#ifndef ERRSTATS_H
#define ERRSTATS_H

#include "pli.h"
#include "contour.h"
#include "cmaes.h"

namespace libcmaes
{
  
  template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class errstats
    {
    public:
    /*- profile likelihood -*/

    /**
     * \brief computes the profile likelihood in dimension k around a previously found optima
     * @param func objective function
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously found optima
     * @param k dimension in which to compute profile likelihood points
     * @param curve whether to store all points during search in order to build a profile likelihood curve
     * @param samplesize number of steps of linesearch in every direction in dimension k
     * @param fup the function deviation for which to compute the profile likelihood
     * @param delta tolerance around fvalue + fup for which to compute the profile likelihood
     * @param maxiters maximum number of linesearch tentatives for computing the profile likelihood
     * @return profile likelihood object
     * @see pli
     */
    static pli profile_likelihood(FitFunc &func,
				  const CMAParameters<TGenoPheno> &parameters,
				  CMASolutions &cmasol,
				  const int &k,
				  const bool &curve=false,
				  const int &samplesize=10,
				  const double &fup=0.1,
				  const double &delta=0.1,
				  const int &maxiters=1e4);

    private:
    /**
     * \brief computes and search the profile likelihood points in dimension k around a previously 
     *        found optima and in a given direction
     * @param func objective function
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously found optima
     * @param k dimension in which to compute profile likelihood points
     * @param neg whether to go on the right (i.e. search direction)
     * @param curve whether to store all points during search in order to build a profile likelihood curve
     * @param samplesize number of steps of linesearch in every direction in dimension k
     * @param fup the function deviation for which to compute the profile likelihood
     * @param delta tolerance around fvalue + fup for which to compute the profile likelihood
     * @param maxiters maximum number of linesearch tentatives for computing the profile likelihood
     */
    static void profile_likelihood_search(FitFunc &func,
					  const CMAParameters<TGenoPheno> &parameters,
					  pli &le,
					  const CMASolutions &cmasol,
					  const int &k,
					  const bool &neg,
					  const int &samplesize,
					  const double &fup,
					  const double &delta,
					  const int &maxiters,
					  const bool &curve);

    /**
     * \brief take a linesearch step in a given direction
     *        Note: the search takes place in geno-space
     * @param func objective function
     * @param parameters stochastic search parameters
     * @param k dimension in which to compute profile likelihood points
     * @param minfvalue current objective function min value
     * @param n number of steps allowed in linesearch
     * @param fup the function deviation for which to compute the profile likelihood
     * @param delta tolerance around fvalue + fup for which to compute the profile likelihood
     * @param linit whether this is the first linesearch call
     * @param eigenve eigenvectors
     * @param d step
     * @param x vector on the line
     */
    static void take_linear_step(FitFunc &func,
				 const CMAParameters<TGenoPheno> &parameters,
				 const int &k,
				 const double &minfvalue,
				 const double &fup,
				 const double &delta,
				 const int &n,
				 const bool &linit,
				 const dMat &eigenve,
				 double &d,
				 dVec &x);

    public:
    /**
     * \brief optimizes an objective function while fixing the value of parameters in several dimensions
     * @param func objective function
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously found optima
     * @param k dimensions in which to fix parameters (i.e. search takes place in all other dimensions)
     * @param vk fixed values of parameters in dimensions of set k
     * @param x0 initial parameter values
     * @param pheno_x0 whether x0 is in phenotype
     * @param pheno_vk whether vk is in phenotype
     * @return optimization solution partial object with a single candidate that is the best candidate in full dimension
     */
    static CMASolutions optimize_vpk(FitFunc &func,
				     const CMAParameters<TGenoPheno> &parameters,
				     const CMASolutions &cmasol,
				     const std::vector<int> &k,
				     const std::vector<double> &vk,
				     const dVec &x0,
				     const bool &pheno_x0=true,
				     const bool &pheno_vk=true);
    
    /**
     * \brief optimizes an objective function while fixing the value of parameters in dimension k
     * @param func objective function
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously found optima
     * @param k dimension into which to fix the parameter (i.e. search takes place in all other dimensions)
     * @param vk fixed value of parameter k
     * @param x0 initial parameter values
     * @param pheno_x0 whether x0 is in phenotype
     * @param pheno_vk whether vk is in phenotype
     * @return optimization solution partial object with a single candidate that is the best candidate in full dimension
     */
    static CMASolutions optimize_pk(FitFunc &func,
				    const CMAParameters<TGenoPheno> &parameters,
				    const CMASolutions &cmasol,
				    const int &k,
				    const double &vk,
				    const dVec &x0,
				    const bool &pheno_x0=true,
				    const bool &pheno_vk=true);
    
    /*- contour -*/
    public:
    /**
     * \brief computes a set of contour points around a function minima, for a deviation fup in objective function value
     * @param func objective function
     * @param px first dimension for contour points computation
     * @param py second dimension for contour points computation
     * @param npoints number of points to be computed in contour
     * @param fup the function deviation for which to compute the contour
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously found optima
     * @param delta tolerance around fvalue + fup for which to compute the profile likelihood
     * @param maxiters maximum number of linesearch tentatives for computing the profile likelihood
     * @return contour object that contains the contour points
     */
    static contour contour_points(FitFunc & func, const int &px, const int &py, const int &npoints, const double &fup,
				  const CMAParameters<TGenoPheno> &parameters,
				  CMASolutions &cmasol,
				  const double &delta=0.1,
				  const int &maxiters=1e4);

    private:
    /**
     * \brief finds crossing point
     * @param parameters stochastic search parameters
     * @param cmasol solution object that contains the previously foundoptima
     * @param fup the function deviation for which to compute the contour
     * @param par pair of dimensions to work in
     * @param pmid middle point in both dimensions
     * @param pdir direction in both dimensions
     * @param ftol tolerance around fvalue + fup
     * @return crossing object
     */
    static fcross cross(FitFunc &func,
			const CMAParameters<TGenoPheno> &parameters,
			CMASolutions &cmasol,
			const double &fup,
			const std::vector<int> &par, const std::vector<double> &pmid,
			const std::vector<double> &pdir, const double &ftol);
    };    
  
}

#endif
