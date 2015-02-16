# CMA-ES, Covariance Matrix Adaptation Evolution Strategy
# Copyright (c) 2015 Inria
# Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
#
# This file is part of libcmaes.
#
# libcmaes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libcmaes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
##

import sys
import numpy as np

# number of static variables at the head of every line (i.e. independent of problem dimension)
single_values = 4

def convert(infile,outprefix):
    # get relevant information from first line
    first_line = ''
    with open(infile,'r') as f:
        first_line = f.readline()
        f.close()
    fspl = first_line.split(' ',3)
    dim = int(fspl[0])
    seed = fspl[1]
    fspl = first_line.split(' / ',1)
    rdate = fspl[1]
    #print 'dim=',dim
    #print 'seed=',seed
    #print 'date=',rdate

    # read data from input file
    dat = np.loadtxt(infile,dtype=float,skiprows=2)   

    # parse data from input file
    fvalue = np.absolute(dat[:,0])
    fevals = dat[:,1]
    sigma = dat[:,2]
    kappa = dat[:,3]
    best_seen_fvalue = dat[:,4]
    median_seen_fvalue = dat[:,5]
    worst_seen_fvalue = dat[:,6]
    min_eigenv = dat[:,7]
    max_eigenv = dat[:,8]
    pos = 9
    best_seen_x = dat[:,range(pos,pos+dim)]
    pos += dim
    eigenv = dat[:,range(pos,pos+dim)]
    pos += dim
    stds = dat[:,range(pos,pos+dim)]
    pos += dim
    xmean = dat[:,range(pos,pos+dim)]
    pos += dim
    
    # create output files and fill them up according to legacy format
    ffit = open(outprefix + 'fit.dat','w')
    ffit.write('% # columns="iteration, evaluation, sigma, axis ratio, bestever, best, median, worst objective function value, further objective values of best, seed=' + seed + ', ' + rdate + '\n')
    fxrecent = open(outprefix + 'xrecentbest.dat','w')
    fxrecent.write('% # iter+eval+sigma+0+fitness+xbest, seed=' + seed + ', ' + rdate + '\n')
    fxmean = open(outprefix + 'xmean.dat','w')
    fxmean.write('% # columns="iteration, evaluation, void, void, void, xmean", seed=' + seed + ', ' + rdate + '\n')
    faxlen = open(outprefix + 'axlen.dat','w')
    faxlen.write('%  columns="iteration, evaluation, sigma, max axis length,  min axis length, all principle axes lengths  (sorted square roots of eigenvalues of C)", seed=' + seed + ', ' + rdate + '\n')
    fstddev = open(outprefix + 'stddev.dat','w')
    fstddev.write('% # columns="iteration, evaluation, sigma, void, void,  stds==sigma*sqrt(diag(C))", seed=' + seed + ', ' + rdate + '\n')

    i = 0
    for fe,si,ka,bs,fv,ms,ws,Me,me in np.nditer([fevals,sigma,kappa,best_seen_fvalue,fvalue,median_seen_fvalue,worst_seen_fvalue,max_eigenv,min_eigenv]):
        ffit.write('%s %s %s %s %s %s %s %s\n' % (i,fe,si,ka,bs,fv,ms,ws))
        xb = ' '.join(map(str,best_seen_x[i:i+1][0]))
        fxrecent.write('%s %s %s %s 0 %s\n' % (i,fe,si,fv,xb))
        xm = ' '.join(map(str,xmean[i:i+1][0]))
        fxmean.write('%s %s 0 0 0 %s\n' % (i,fe,xm))
        ei = ' '.join(map(str,eigenv[i:i+1][0]))
        faxlen.write('%s %s %s %s %s %s\n' % (i,fe,si,Me,me,ei))
        sd = ' '.join(map(str,stds[i:i+1][0]))
        fstddev.write('%s %s %s 0 0 %s\n' % (i,fe,si,sd))
        i = i + 1
    ffit.close()
    fxrecent.close()
    fxmean.close()
    faxlen.close()
    fstddev.close()

if __name__ == "__main__":
    convert(sys.argv[1],'outcmaes')
