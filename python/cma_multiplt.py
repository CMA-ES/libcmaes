#!/usr/bin/env python
"""In a OS shell::

    python cma_multiplt.py data_file_name
    
or in a python shell::

    import cma_multiplt as lcmaplt
    lcmaplt.plot(data_file_name)
    
"""
# CMA-ES, Covariance Matrix Adaptation Evolution Strategy
# Copyright (c) 2014 Inria
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

import sys, pylab, csv
import numpy as np
from matplotlib.pylab import subplot, semilogy, grid, title
# from matplotlib.pylab import figure, subplot, semilogy, hold, grid, axis, title, text, xlabel, isinteractive, draw, gcf
# TODO: the above direct imports clutter the interface in a Python shell

# number of static variables at the head of every line (i.e. independent of problem dimension)
single_values = 4

def plot(filename):
    # read data into numpy array
    dat = np.loadtxt(filename,dtype=float)

    dim = int(np.ceil(np.shape(dat)[1] - single_values) / 3) # we estimate the problem dimension from the data
    #print dim

    fvalue = np.absolute(dat[:,0])
    fevals = dat[:,1]
    sigma = dat[:,2]
    kappa = dat[:,3]
    if dim > 0:
        eigenvc = []
        for c in range(single_values,single_values+dim):
            eigenvc.append(c)
        eigenv = dat[:,eigenvc]
        stdsc = []
        for c in range(single_values+dim,single_values+2*dim):
            stdsc.append(c)
        stds = dat[:,stdsc]
        minstds = np.amin(stds,axis=1)
        maxstds = np.amax(stds,axis=1)
        xmeanc = []
        for c in range(single_values+2*dim,single_values+3*dim):
            xmeanc.append(c)
        xmean = dat[:,xmeanc]

    # plot data.
    pylab.rcParams['font.size'] = 10
    xlab = "function evaluations"

    # plot fvalue, sigma, kappa
    if dim > 0:
        subplot(221)
    semilogy(fevals,fvalue,'b')
    semilogy(fevals,sigma,'g')
    semilogy(fevals,kappa,'r')
    if dim > 0:
        semilogy(fevals,sigma*minstds,'y')
        semilogy(fevals,sigma*maxstds,'y')
    title('f-value (blue), sigma (green), kappa (red)')
    grid(True)

    if dim == 0:
        pylab.xlabel(xlab)
        pylab.show();
        msg = '  --- press return to continue --- '
        raw_input(msg) if sys.version < '3' else input(msg)
        sys.exit(1)

    # plot xmean
    subplot(222)
    pylab.plot(fevals,xmean)
    title('Object Variables (mean, ' + str(dim) + '-D)')
    grid(True)

    # plot eigenvalues
    subplot(223)
    semilogy(fevals,eigenv,'-b')
    pylab.xlabel(xlab)
    title('Eigenvalues')
    grid(True)

    # plot std deviations
    subplot(224)
    semilogy(fevals,stds)
    pylab.xlabel(xlab)
    title('Standard Deviation in all coordinates')
    grid(True)

    pylab.ion()
    pylab.draw()
    pylab.show()

if __name__ == "__main__":
    plot(sys.argv[1])
    msg = '  --- press return to continue --- '
    raw_input(msg) if sys.version < '3' else input(msg)
