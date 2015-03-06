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

#include "errstats.h"
#include <llogging.h>
#include <iostream>

namespace libcmaes
{
  
  template <class TGenoPheno>
  pli errstats<TGenoPheno>::profile_likelihood(FitFunc &func,
					       const CMAParameters<TGenoPheno> &parameters,
					       CMASolutions &cmasol,
					       const int &k,
					       const bool &curve,
					       const int &samplesize,
					       const double &fup,
					       const double &delta,
					       const int &maxiters)
  {
    dVec x = cmasol.best_candidate().get_x_dvec();
    double minfvalue = cmasol.best_candidate().get_fvalue();
    dVec phenox = parameters._gp.pheno(x);
    
    //debug
    /*std::cout << "xk=" << x[k] << " / minfvalue=" << minfvalue << std::endl;
      std::cout << "phenox=" << phenox.transpose() << std::endl;*/
    //debug

    pli le(k,samplesize,parameters._dim,parameters._gp.pheno(x),minfvalue,fup,delta);
    
    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,false,samplesize,fup,delta,maxiters,curve); // positive direction
    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,true,samplesize,fup,delta,maxiters,curve);  // negative direction
    
    le.setErrMinMax();
    cmasol._pls.insert(std::pair<int,pli>(k,le));
    return le;
  }

  template <class TGenoPheno>
  void errstats<TGenoPheno>::profile_likelihood_search(FitFunc &func,
						       const CMAParameters<TGenoPheno> &parameters,
						       pli &le,
						       const CMASolutions &cmasol,
						       const int &k,
						       const bool &neg,
						       const int &samplesize,
						       const double &fup,
						       const double &delta,
						       const int &maxiters,
						       const bool &curve)
  {
    int sign = neg ? -1 : 1;
    dVec x = cmasol.best_candidate().get_x_dvec();
    double xk = x[k];
    double minfvalue = cmasol.best_candidate().get_fvalue();
    double nminfvalue = minfvalue;
    CMASolutions citsol = cmasol;
    int i = 0;
    int n = 10;
    double d = sign * xk * 0.1; // adhoc.
    bool linit = true;
    while(true)
      {
	// get a new xk point.
	errstats<TGenoPheno>::take_linear_step(func,parameters,k,minfvalue,fup,delta,n,linit,cmasol.eigenvectors(),d,x);
	linit = false;
	
	//debug
	//std::cout << "new xk point: " << x.transpose() << std::endl;
	//debug
	
	// minimize.
	CMASolutions ncitsol = errstats<TGenoPheno>::optimize_pk(func,parameters,citsol,k,x[k],x,false,false);
	if (ncitsol._run_status < 0)
	  {
	    LOG(ERROR) << "profile likelihood linesearch: optimization error " << ncitsol._run_status << " -- " << ncitsol.status_msg() << std::endl;
	    // pad and return.
	    for (int j=0;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = le._fvaluem[samplesize];
		le._xm.row(samplesize+sign*(1+j)) = le._xm.row(samplesize);
		le._err[samplesize+sign*(1+j)] = ncitsol._run_status;
	      }
	      return;
	  }
	else if (ncitsol.best_candidate().get_fvalue() < 1e-1*minfvalue)
	  {
	    LOG(ERROR) << "profile likelihood finds new minimum: " << ncitsol.best_candidate().get_fvalue() << std::endl;
	    // pad and return.
	    for (int j=0;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = le._fvaluem[samplesize];
		le._xm.row(samplesize+sign*(1+j)) = le._xm.row(samplesize);
		le._err[samplesize+sign*(1+j)] = ncitsol._run_status;
	      }
	      return;
	  }
	else // update current point and solution.
	  {
	    x = ncitsol.best_candidate().get_x_dvec();
	    nminfvalue = ncitsol.best_candidate().get_fvalue();
	  }
	
	// store points.
	dVec phenobx = parameters._gp.pheno(ncitsol.best_candidate().get_x_dvec());
	if (curve)
	  {
	    le._fvaluem[samplesize+sign*(1+i)] = ncitsol.best_candidate().get_fvalue();
	    le._xm.row(samplesize+sign*(1+i)) = phenobx.transpose();
	    le._err[samplesize+sign*(1+i)] = ncitsol._run_status;
	  }

	bool iterend = (fabs(nminfvalue-minfvalue-fup) <= 0.1 * fup);

	//debug
	//std::cerr << "iterend=" << iterend << " / nminfvalue=" << nminfvalue << std::endl;
	//debug
	
	if (!curve && iterend)
	  {
	    // pad and return.
	    for (int j=0;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = ncitsol.best_candidate().get_fvalue();
		le._xm.row(samplesize+sign*(1+j)) = phenobx.transpose();
		le._err[samplesize+sign*(1+j)] = ncitsol._run_status;
	      }
	    return;
	  }
	++i;
	if (curve && i == samplesize)
	  break;
	if (i == maxiters)
	  break;
      }
  }

  template <class TGenoPheno>
  void errstats<TGenoPheno>::take_linear_step(FitFunc &func,
					      const CMAParameters<TGenoPheno> &parameters,
					      const int &k,
					      const double &minfvalue,
					      const double &fup,
					      const double &delta,
					      const int &n,
					      const bool &linit,
					      const dMat &eigenve,
					      double &d,
					      dVec &x)
  {
    double fdelta = delta * fup;
    dVec xtmp = x;
    dVec phenoxtmp = parameters._gp.pheno(xtmp);
    double fvalue = func(phenoxtmp.data(),xtmp.size());
    double fdiff = fabs(fvalue - minfvalue);
    
    //debug
    //std::cerr << "d=" << d << " / fdiff=" << fdiff << " / xtmpk=" << xtmp[k] << std::endl;
    //debug
    
    dVec dv = dVec::Zero(parameters.dim());
    dv[k] = d;
    int i = n;
    if (fdiff > fup + fdelta) // above
      {
	x -= eigenve.transpose() * dv;
	while (fdiff > fup + fdelta
	       && i > 0
	       && phenoxtmp[k] >= parameters._gp.get_boundstrategy().getPhenoLBound(k)
	       && phenoxtmp[k] <= parameters._gp.get_boundstrategy().getPhenoUBound(k))
	  {
	    d /= 2.0;
	    dv[k] = d;
	    xtmp = x + eigenve.transpose() * dv; // search in rotated space
	    phenoxtmp = parameters._gp.pheno(xtmp);
	    fvalue = func(phenoxtmp.data(),xtmp.size());
	    fdiff = fabs(fvalue - minfvalue);

	    //debug
	    //std::cerr << "xtmpk=" << xtmp[k] << " / fvalue=" << fvalue << " / fdiff=" << fdiff << " / dxk=" << d << std::endl;
	    //debug

	    if (!linit)
	      --i;
	  }
      }
    else // below
      {
	while (fdiff < fup - fdelta
	       && i > 0
	       && phenoxtmp[k] >= parameters._gp.get_boundstrategy().getPhenoLBound(k)
	       && phenoxtmp[k] <= parameters._gp.get_boundstrategy().getPhenoUBound(k))
	  {
	    d *= 2.0;
	    dv[k] = d;
	    xtmp = x + eigenve.transpose() * dv; // search in rotated space
	    phenoxtmp = parameters._gp.pheno(xtmp);
	    fvalue = func(phenoxtmp.data(),xtmp.size());
	    fdiff = fabs(fvalue - minfvalue);

	    //debug
	    //std::cerr << "xtmpk=" << xtmp[k] << " / fvalue=" << fvalue << " / fdiff=" << fdiff << " / dxk=" << d << std::endl;
	    //debug

	    if (!linit) // first pass must reach above threshold.
	      --i;
	  }
      }
    x[k] = xtmp[k];
  }
  
  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_vpk(FitFunc &func,
						  const CMAParameters<TGenoPheno> &parameters,
						  const CMASolutions &cmasol,
						  const std::vector<int> &k,
						  const std::vector<double> &vk,
						  const dVec &x0,
						  const bool &pheno_x0,
						  const bool &pheno_vk)
  {
    dVec rx0;
    if (!pheno_x0)
      rx0 = parameters.get_gp().pheno(x0);
    else rx0 = x0;
    double sigma = 0.0;
    TGenoPheno ngp = parameters.get_gp();
    ngp.remove_dimensions(k);
    for (int i=k.size()-1;i>=0;i--)
      {
	sigma = std::max(sigma,fabs(cmasol._candidates.at(0).get_x_dvec()[k[i]]-vk[i])); // heuristic: sigma as max distance to minima
	removeElement(rx0,k[i]);
      }
    
    // transform the point of interest from genotype to phenotype
    dVec pvk = dVec::Zero(x0.size());
    for (size_t i=0;i<vk.size();i++)
      pvk(k[i]) = vk.at(i);
    if (!pheno_vk)
      pvk = parameters.get_gp().pheno(pvk);
     
    if (rx0.size() == 0) // if all variables are fixed, simply return value for this point
      {
	double fval = func(pvk.data(),pvk.size());
	CMASolutions rcmasol;
	rcmasol._candidates.emplace_back(fval,pvk);
	rcmasol._best_candidates_hist.push_back(cmasol._candidates.at(0));
	rcmasol._nevals++;
	return rcmasol;
      }
    CMAParameters<TGenoPheno> nparameters(rx0.size(),rx0.data(),sigma,parameters.lambda(),
					  parameters.get_seed(),ngp);
    nparameters.set_initial_fvalue(true);
    nparameters.set_ftarget(parameters.get_ftarget());
    //nparameters.set_quiet(false);
    
    FitFunc rfunc = [func,k,pvk](const double *x, const int N)
      {
	dVec nx(N);
	for (int i=0;i<N;i++)
	  nx[i] = x[i];
	for (size_t i=0;i<k.size();i++)
	  {
	    addElement(nx,k[i],pvk[k[i]]); // in phenotype
	  }
	return func(nx.data(),nx.size());
      };
        
    CMASolutions cms = cmaes<TGenoPheno>(rfunc,nparameters);
    dVec nx = cms.best_candidate().get_x_dvec();
    if (!pheno_vk)
      {
	for (size_t i=0;i<k.size();i++)
	  addElement(nx,k[i],vk[i]);
      }
      else
      {
	// XXX: wish to not have to rely on the geno() function, which 
	//      remains unused elsewhere.
	//      profile likelihood proceeds in genotype space, whereas contour
	//      proceeds in phenotype space, therefore the present function 
	//      is required to work in both cases.
	dVec gpvk = parameters.get_gp().geno(pvk);
	for (size_t i=0;i<k.size();i++)
	  addElement(nx,k[i],gpvk[k[i]]);
      }
    CMASolutions rcmasol;
    rcmasol._candidates.emplace_back(cms.best_candidate().get_fvalue(),nx); // in genotype
    rcmasol._best_candidates_hist.push_back(rcmasol._candidates.at(0));
    rcmasol._nevals = cms._nevals;
    rcmasol._run_status = cms.run_status();
    return rcmasol;
  }

  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_pk(FitFunc &func,
						 const CMAParameters<TGenoPheno> &parameters,
						 const CMASolutions &cmasol,
						 const int &k,
						 const double &vk,
						 const dVec &x0,
						 const bool &pheno_x0,
						 const bool &pheno_vk)
  {
    std::vector<int> tk = {k};
    std::vector<double> tvk = {vk};
    return errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,tk,tvk,x0,pheno_x0,pheno_vk);
  }

  template <class TGenoPheno>
  contour errstats<TGenoPheno>::contour_points(FitFunc & func, const int &px, const int &py, const int &npoints, const double &fup,
					       const CMAParameters<TGenoPheno> &parameters,
					       CMASolutions &cmasol,
					       const double &delta,
					       const int &maxiters)
  {
    // find first two points.
    int samplesize = 10;
    pli plx,ply;
    if (!cmasol.get_pli(px,plx))
      {
	errstats<TGenoPheno>::profile_likelihood(func,parameters,cmasol,px,false,samplesize,fup,delta,maxiters);
	cmasol.get_pli(px,plx); // in phenotype
      }
    
    // find second two points.
    if (!cmasol.get_pli(py,ply))
      {
	errstats<TGenoPheno>::profile_likelihood(func,parameters,cmasol,py,false,samplesize,fup,delta,maxiters);
	cmasol.get_pli(py,ply); // in phenotype
      }

    dVec phenox = cmasol.best_candidate().get_x_pheno_dvec(parameters); // in phenotype
    double valx = phenox(px);
    double valy = phenox(py);
    
    // find upper y value for x parameter.
    CMASolutions exy_up = errstats<TGenoPheno>::optimize_pk(func,parameters,cmasol,px,valx+plx._errmax,parameters.get_x0min());
    //std::cout << "exy_up=" << exy_up.best_candidate().get_x_dvec().transpose() << std::endl;
    
    // find lower y value for x parameter.
    CMASolutions exy_lo = errstats<TGenoPheno>::optimize_pk(func,parameters,cmasol,px,valx+plx._errmin,parameters.get_x0min());
    //std::cout << "exy_lo=" << exy_lo.best_candidate().get_x_dvec().transpose() << std::endl;
    
    // find upper x value for y parameter.
    CMASolutions eyx_up = errstats<TGenoPheno>::optimize_pk(func,parameters,cmasol,py,valy+ply._errmax,parameters.get_x0min());
    //std::cout << "eyx_up=" << eyx_up.best_candidate().get_x_dvec().transpose() << std::endl;
    
    // find lower x value for y parameter.
    CMASolutions eyx_lo = errstats<TGenoPheno>::optimize_pk(func,parameters,cmasol,py,valy+ply._errmin,parameters.get_x0min());
    //std::cout << "eyx_lo=" << eyx_lo.best_candidate().get_x_dvec().transpose() << std::endl;
    
    // early contour in phenotype
    contour c;
    c.add_point(valx+plx._errmin,exy_lo.best_candidate().get_x_pheno_dvec(parameters)(py));
    c.add_point(eyx_lo.best_candidate().get_x_pheno_dvec(parameters)(px),valy+ply._errmin); 
    c.add_point(valx+plx._errmax,exy_up.best_candidate().get_x_pheno_dvec(parameters)(py));
    c.add_point(eyx_up.best_candidate().get_x_pheno_dvec(parameters)(px),valy+ply._errmax);

    double scalx = 1.0/(plx._errmax - plx._errmin); // in phenotype
    double scaly = 1.0/(ply._errmax - ply._errmin);

    //debug
    //std::cout << "contour:" << c << std::endl;
    //debug
    
    // more than 4 points.
    for (int i=4;i<npoints;i++)
      {
	//debug
	//std::cout << "=> generating point #" << i << std::endl;
	//debug
	
	//- check on max budget, and return if exceeded.
	
	//- get most distant points, in phenotype.
	std::vector<std::pair<double,double>>::iterator idist1 = c._points.end()-1;
	std::vector<std::pair<double,double>>::iterator idist2 = c._points.begin();
	double dx = idist1->first - idist2->first;
	double dy = idist1->second - idist2->second;
	double bigdis = scalx*scalx*dx*dx + scaly*scaly*dy*dy;

	for (auto ipair = c._points.begin();ipair!=c._points.end()-1;ipair++)
	  {
	    dx = ipair->first - (ipair+1)->first;
	    dy = ipair->second - (ipair+1)->second;
	    double dist = scalx*scalx*dx*dx + scaly*scaly*dy*dy;
	    if (dist > bigdis)
	      {
		bigdis = dist;
		idist1 = ipair;
		idist2 = ipair+1;
	      }
	    //std::cout << "idist10=" << idist1->first << " -- idist20=" << idist2->first << std::endl;
	  }
	
	// - select mid-range point x and direction dir along the two axis of interest.
	double a1 = 0.5;
	double a2 = 0.5;
	double sca = 1.0;
	double xmidcr = a1*idist1->first + a2*idist2->first;
	double ymidcr = a1*idist1->second + a2*idist2->second;
	double xdir = idist2->second - idist1->second;
	double ydir = idist1->first - idist2->first;
	double scalfac = sca*std::max(fabs(xdir*scalx), fabs(ydir*scaly));
	double xdircr = xdir/scalfac;
	double ydircr = ydir/scalfac;
	std::vector<double> pmid = {xmidcr, ymidcr};
	std::vector<double> pdir = {xdircr, ydircr};
	std::vector<int> par = {px,py};
	
	//debug
	/*std::cout << "idist10=" << idist1->first << " / idist11=" << idist1->second << " / idist20=" << idist2->first << " / idist21=" << idist2->second << std::endl;
	  std::cout << "pmid0=" << pmid[0] << " / pmid1=" << pmid[1] << " / pdir0=" << pdir[0] << " / pdir1=" << pdir[1] << std::endl;*/
	//debug
	
	// find crossing point from x with direction dir where function is equal to min + fup.
	fcross fc = errstats<TGenoPheno>::cross(func,parameters,cmasol,fup,par,pmid,pdir,parameters._ftolerance);
	if (fc._nevals == 0.0) // dummy
	  continue;
	
	//debug
	//std::cout << "fcross point=" << fc._x.transpose() << std::endl;
	//debug

	if (idist2 == c._points.begin())
	  {
	    c.add_point(fc._x(par[0]),fc._x(par[1]));
	  }
	else
	  {
	    c.add_point(idist2,fc._x(par[0]),fc._x(par[1]));
	  }
      }

    //debug
    //std::cout << "number of contour points=" << c._points.size() << std::endl;
    //debug
    
    return c;
  }

  template <class TGenoPheno>
  fcross errstats<TGenoPheno>::cross(FitFunc &func,
				     const CMAParameters<TGenoPheno> &parameters,
				     CMASolutions &cmasol,
				     const double &fup,
				     const std::vector<int> &par, const std::vector<double> &pmid,
				     const std::vector<double> &pdir, const double &ftol)
  {
    double aopt = 0.0;
    std::vector<double> alsb(3,0.0), flsb(3,0.0);
    std::vector<CMASolutions> cmasols;
    double aminsv = cmasol.best_candidate().get_fvalue();
    double aim = aminsv + fup;
    double tla = ftol;
    double tlf = ftol*fup;
    int nfcn = 0;
    unsigned int maxitr = 15, ipt = 0;
    dVec nx;

    //debug
    //std::cout << "aminsv=" << aminsv << " / aim=" << aim << " / tla=" << tla << " / tlf=" << tlf << std::endl;
    //debug
    
    // get a first optimized point.
    CMASolutions cmasol1 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmid,parameters.get_x0min(),true,true); //TODO: result has wrong best candidate x set as phenotype...
    alsb[0] = 0.0;
    flsb[0] = cmasol1.best_candidate().get_fvalue();
    flsb[0] = std::max(flsb[0],aminsv+0.1*fup);
    nfcn += cmasol1._nevals;
    cmasols.push_back(cmasol1);
    ipt++;

    //debug
    //std::cout << "contour / fvalue=" << cmasol1.best_candidate().get_fvalue() << " / optimized point1=" << cmasol1.best_candidate().get_x_pheno_dvec(parameters).transpose() << std::endl;
    //debug
    
    // update aopt and get a second optimized point.
    aopt = sqrt(fup/(flsb[0]-aminsv))-1.0;
    if (aopt > 1.0)
      aopt = 1.0;
    else if (aopt < -0.5)
      aopt = 0.5;
    std::vector<double> pmiddir2 = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
    CMASolutions cmasol2 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmiddir2,parameters.get_x0min());
    alsb[1] = aopt;
    flsb[1] = cmasol2.best_candidate().get_fvalue();
    nfcn += cmasol2._nevals;
    double dfda = (flsb[1]-flsb[0])/(alsb[1]-alsb[0]);
    cmasols.push_back(cmasol2);
    ipt++;

    //debug
    //std::cout << "contour / fvalue=" << cmasol2.best_candidate().get_fvalue() << " / optimized point2=" << cmasol2.best_candidate().get_x_dvec().transpose() << std::endl;
    //debug
    
    // study slope between the two points.
  slope:
    if (dfda < 0.0)
      {

	//debug
	//std::cout << "negative slope\n";
	//debug
	
	// if negative slope, update until positive.
	unsigned int maxlk = maxitr - ipt;
	for (unsigned int it=0;it<maxlk;it++)
	  {
	    alsb[0] = alsb[1];
	    flsb[0] = flsb[1];
	    aopt = alsb[0] + 0.2*(it+1);
	    std::vector<double> pmidt = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
	    CMASolutions cmasolt = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmidt,parameters.get_x0min());
	    alsb[1] = aopt;
	    flsb[1] = cmasolt.best_candidate().get_fvalue();
	    dfda = (flsb[1]-flsb[0])/(alsb[1]-alsb[0]);
	    nfcn += cmasolt._nevals;
	    cmasols.at(1) = cmasolt;
	    if (dfda > 0.0)
	      break;
	  }
      }

    // once positive slope, find a third point.
    int tentatives = 0;
  point3:
    //debug
    //std::cout << "positive slope, dfda=" << dfda << std::endl;
    //debug
    
    tentatives++;
    aopt = alsb[1] + (aim-flsb[1])/dfda;
    double fdist = std::min(fabs(aim  - flsb[0]), fabs(aim  - flsb[1]));
    double adist = std::min(fabs(aopt - alsb[0]), fabs(aopt - alsb[1]));
    if (fabs(aopt) > 1.0)
      tla = ftol*fabs(aopt);
    if (adist < tla && fdist < tlf)
      {
	// return second optimized point.
	//debug
	//std::cout << "below tolerance, returning second optimized point\n";
	//debug
	return fcross(cmasol2.best_candidate().get_fvalue(),
		      nfcn,cmasol2.best_candidate().get_x_pheno_dvec(parameters));
      }
    double bmin = std::min(alsb[0], alsb[1]) - 1.;
    if (aopt < bmin) aopt = bmin;
    double bmax = std::max(alsb[0], alsb[1]) + 1.;
    if (aopt > bmax) aopt = bmax;

    // get third point.
    std::vector<double> pmiddir3 = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
    CMASolutions cmasol3 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmiddir3,parameters.get_x0min());
    alsb[2] = aopt;
    flsb[2] = cmasol3.best_candidate().get_fvalue();
    nfcn += cmasol3._nevals;
    if (cmasols.size() < 3)
      {
	cmasols.push_back(cmasol3);
      }
    else 
      {
	cmasols.at(2) = cmasol3;
      }
    
    //debug
    //std::cout << "contour / fvalue=" << cmasol3.best_candidate().get_fvalue() << " / optimized point3=" << cmasol3.best_candidate().get_x_dvec().transpose() << std::endl;
    //debug
    
    // from three points < or > objective, decide what to do.
    double ecarmn = fabs(flsb[2] - aim);
    double ecarmx = 0.;
    unsigned int ibest = 2;
    unsigned int iworst = 0;
    unsigned int noless = 0;

    for(unsigned int i = 0; i < 3; i++)
      {
	double ecart = fabs(flsb[i] - aim);
	if(ecart > ecarmx) {
	  ecarmx = ecart;
	  iworst = i;
	}
	if(ecart < ecarmn)
	  {
	    ecarmn = ecart;
	    ibest = i;
	  }
	if(flsb[i] < aim) noless++;
      }
    
    // at least one on each side of AIM (contour)
    if(noless == 1 || noless == 2)
      {
	// XXX: could do parabola instead.
	int srefp = -1;
	bool bestbelow = (flsb[ibest] < aim);
	double srefecar = ecarmx;
	for (int i=0;i<3;i++)
	  {
	    if (i != (int)ibest)
	      {
		if ((bestbelow && flsb[i] > aim)
		    || (!bestbelow && flsb[i] < aim))
		  {
		    double ecar = fabs(flsb[i]-aim);
		    if (ecar <= srefecar)
		      {
			srefecar = ecar;
			srefp = i;
		      }
		  }
	      }
	  }
	dVec srefp_phenox = cmasols.at(srefp).best_candidate().get_x_pheno_dvec(parameters);
	dVec ibest_phenox = cmasols.at(ibest).best_candidate().get_x_pheno_dvec(parameters);
	int sdir0 = (srefp_phenox(par[0])-ibest_phenox(par[0]) > 0.0) ? 1 : -1;
	int sdir1 = (srefp_phenox(par[1])-ibest_phenox(par[1]) > 0.0) ? 1 : -1;
	dVec xstart;
	double nminfvalue;
	if (cmasols.at(ibest).best_candidate().get_fvalue() < aim)
	  {
	    xstart = ibest_phenox;
	    nminfvalue = cmasols.at(ibest).best_candidate().get_fvalue();
	  }
	else
	  {
	    xstart = srefp_phenox;
	    nminfvalue = cmasols.at(srefp).best_candidate().get_fvalue();
	  }
	double dxk0 = sdir0 * fabs(xstart(par[0])) * 0.1;
	double dxk1 = sdir1 * fabs(xstart(par[1])) * 0.1;

	//debug
	//std::cout << "sdir0=" << sdir0 << " / sdir1=" << sdir1 << " / dxk0=" << dxk0 << " / dxk1=" << dxk1 << std::endl;
	//debug
	
	CMASolutions citsol = cmasols.at(ibest);
	while (true)
	  {
	    // advance incrementally
	    xstart(par[0]) += dxk0;
	    xstart(par[1]) += dxk1;
	    
	    //debug
	    //std::cout << "xstart=" << xstart.transpose() << std::endl;
	    //debug
	    
	    std::vector<double> vxk = {xstart(par[0]),xstart(par[1])};
	    CMASolutions ncitsol = errstats<TGenoPheno>::optimize_vpk(func,parameters,citsol,par,vxk,parameters.get_x0min());
	    if (ncitsol._run_status < 0)
	      {
		LOG(WARNING) << "contour linesearch: optimization error " << ncitsol._run_status << std::endl;
	      }
	    nminfvalue = ncitsol.best_candidate().get_fvalue();
	    if (nminfvalue > aim)
	      break;
	    citsol = ncitsol;
	  }

	//debug
	//std::cout << "contour linesearch best point=" << citsol.best_candidate().get_x_dvec().transpose() << std::endl;
	//debug

	return fcross(citsol._candidates.at(0).get_fvalue(),nfcn,citsol._candidates.at(0).get_x_pheno_dvec(parameters));
      }
    // if all three points are above aim, third point must be the closest to AIM, return it
    if(noless == 0 && ibest != 2)
      {
	//debug
	//std::cout << "all points above, returning closest point as best\n";
	//debug

	return fcross(flsb[ibest],
		      nfcn,cmasols[ibest].best_candidate().get_x_pheno_dvec(parameters));
      }
    // if all three below and third is not best then the slope has again gone negative,
    // re-iterate and look for positive slope
    if(noless == 3 && ibest != 2)
      {
	//debug
	//std::cout << "slope is again negative\n";
	//debug
	
	alsb[1] = alsb[2];
	flsb[1] = flsb[2];
	goto slope;
      }

    // in other case new straight line thru first two points
    //debug
    //std::cout << "new straight line through first two points\n";
    //debug

    // restart from point 3...
    flsb[iworst] = flsb[2];
    alsb[iworst] = alsb[2];
    dfda = (flsb[1] - flsb[0])/(alsb[1] - alsb[0]);
    if (tentatives < 10)
      goto point3;

    LOG(WARNING) << "returning empty cross point when drawing contour\n";
    return fcross();
  }
  
  template class errstats<GenoPheno<NoBoundStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy>>;
  template class errstats<GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
