#!/usr/bin/env python

import sys
import numpy as np
import pylab as py
from scipy.interpolate import interp1d
from numpy.polynomial.chebyshev import *

_, base_function = sys.argv

N, xi_cut = int(1e3), 6.
xi_points = np.linspace(0, xi_cut, N+1)
xi = 0.5*(xi_points[0:N] + xi_points[1:N+1])

eta, A_, B_, _, _, _ = np.loadtxt("./tables/ci-functions.txt").T
eta, A_ = np.loadtxt("./tables/A.txt").T
eta, B_ = np.loadtxt("./tables/B.txt").T

data = {
    'A': A_, 'B': B_
}[base_function]
base = interp1d(eta, data, kind='cubic')

def plot(func, label):
    py.plot(xi, func(xi)-base(xi), label=label)

#plot(interp1d(eta, data, kind='linear'), 'linear')
#plot(interp1d(eta, data, kind='quadratic'), 'quadratic')
#plot(np.poly1d(np.polyfit(eta, data, 16)), 'polyfit16')
#plot(Chebyshev(chebfit(eta, data, 16)), 'chebfit16')
#plot(lambda x: np.poly1d(np.polyfit(eta**2, data, 20))(x**2), 'quad_polyfit10')
for i in xrange(18, 43):
    if not i % 5:
        plot(np.poly1d(np.polyfit(eta, data, i)), 'polyfit{0}'.format(i))
#    plot(lambda x: Chebyshev(chebfit(eta**2, data, i, w=eta**.5))(x**2), 'quad_chebfit{0}'.format(i))

py.axhline(0, c='k', ls=':', lw=.5)
for x in eta:
    py.axvline(x, c='k', ls=':', lw=.5)
py.legend()
py.show()

