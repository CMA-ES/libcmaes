##
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
from matplotlib.pylab import figure, ioff, ion, subplot, semilogy, hold, plot, grid, axis, title, text, xlabel, isinteractive, draw, gcf
from numpy import *

# number of static variables at the head of every line (i.e. independent of problem dimension)
single_values = 4

# read data into numpy array
dat = loadtxt(sys.argv[1],dtype=float)

dim = int(ceil(np.shape(dat)[1] - single_values - 3) / 3) # we estimate the problem dimension from the data
#print dim

fvalue = np.absolute(dat[:,0])
fevals = dat[:,1]
sigma = dat[:,2]
kappa = dat[:,3]
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
trainerr = dat[:,single_values+3*dim+1]
testerr = dat[:,single_values+3*dim+2]
stesterr = dat[:,single_values+3*dim+3]

# plot data.
pylab.rcParams['font.size'] = 10
xlab = "function evaluations"

# plot fvalue, sigma, kappa
subplot(321)
semilogy(fevals,fvalue,'b')
semilogy(fevals,sigma,'g')
semilogy(fevals,kappa,'r')
semilogy(fevals,sigma*minstds,'y')
semilogy(fevals,sigma*maxstds,'y')
title('f-value (blue), sigma (green), kappa (red)')
grid(True)

# plot xmean
subplot(322)
plot(fevals,xmean)
title('Object Variables (mean, ' + str(dim) + '-D)')
grid(True)

# plot eigenvalues
subplot(323)
semilogy(fevals,eigenv,'-b')
pylab.xlabel(xlab)
title('Eigenvalues')
grid(True)

# plot std deviations
subplot(324)
semilogy(fevals,stds)
pylab.xlabel(xlab)
title('Standard Deviation in all coordinates')
grid(True)

# plot std deviations
subplot(325)
pylab.ylim(0,1)
plot(fevals,trainerr,'r')
plot(fevals,testerr,'g')
plot(fevals,stesterr,'b')
title('Surrogate error: train (red), test (green), smoothed test (blue)')
grid(True)

pylab.show()

msg = '  --- press return to continue --- '
raw_input(msg) if sys.version < '3' else input(msg)
