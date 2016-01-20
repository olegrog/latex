#!/usr/bin/env python

import sys
import pylab as py
import numpy as np
import numpy.random
from functools import partial

korobov3 = [
    [ 10007,    [ 1, 544,     5733    ]],
    [ 50021,    [ 1, 12962,   42926   ]],
    [ 100003,   [ 1, 47283,   15021   ]],
    [ 500009,   [ 1, 33606,   342914  ]],
    [ 1000003,  [ 1, 342972,  439897  ]],
    [ 2000003,  [ 1, 235672,  1208274 ]],
    [ 10000019, [ 1, 4341869, 594760  ]]
]

korobov5 = [
    [ 100003,   [ 1, 11729,   65316,   68384,   51876   ]],
    [ 200003,   [ 1, 62638,   60193,   112581,  142904  ]],
    [ 500009,   [ 1, 191775,  488648,  283447,  69999   ]],
    [ 1000003,  [ 1, 335440,  656043,  403734,  126676  ]],
    [ 2000003,  [ 1, 701679,  680513,  965077,  1248525 ]],
    [ 10000019, [ 1, 3669402, 5455092, 7462912, 2188321 ]],
]

aver = int(sys.argv[1]) if len(sys.argv) > 1 else 1
dim = 5
find_N = lambda N: [line[0] for line in korobov5 if line[0] >= N][0]
baker = lambda x: 1 - np.fabs(2*x-1)

def get_grid1(N):
    return np.random.rand(N, dim).T

def get_grid2(N, transform=lambda x: x, shuffle=False):
    N, coeffs = [line for line in korobov5 if line[0] == N][0]
    grid = transform(np.modf(np.outer(1.+np.arange(N), coeffs) / N + np.random.rand(dim))[0])
    if shuffle:
        np.random.shuffle(grid)
    return grid.T
    #return grid[np.sum(np.square(2*grid[:,0:3]-1), axis=1) < 1].T

#func1 = lambda x: np.cos(.5*np.pi*x)*np.pi/2
cut = 5.4
maxw = lambda x: np.exp(-(cut*(2*x-1))**2)*2*cut/np.pi**.5
maxw_r = lambda x: x**2*np.exp(-(cut*x)**2)*4*cut**3/np.pi**.5
maxw_r6 = lambda x: (2*x-1)**6*np.exp(-(cut*(2*x-1))**2)*16*cut**7/np.pi**.5/15
sin_pi = lambda x: np.sin(np.pi*x)/2*np.pi
func5_car = lambda x: maxw(x[0])*maxw(x[1])*maxw(x[2])*sin_pi(x[3])
func5_sph = lambda x: maxw_r(x[0])*sin_pi(x[1])*sin_pi(x[3])

X = np.array(map(lambda x: x[0], korobov5), dtype=np.float128)
MC_sph, MC_car, QMC_sph, QMC_car = np.copy(X), np.copy(X), np.copy(X), np.copy(X)

calc_value = lambda func, grid, N: np.sum(func(grid(N)))/N
average = lambda func, aver: np.average([ func() for i in xrange(aver)])
for i, N in enumerate(X):
    MC_sph[i] = average(partial(calc_value, func5_sph, get_grid1, N), aver)
    MC_car[i] = average(partial(calc_value, func5_car, get_grid1, N), aver)
    QMC_sph[i] = average(partial(calc_value, func5_sph, get_grid2, N), aver)
    QMC_car[i] = average(partial(calc_value, func5_car, get_grid2, N), aver)

err = lambda x: np.fabs(x - 1)
c = aver**-.5 
py.plot(X, err(MC_sph), label='MC Spherical')
py.plot(X, err(MC_car), label='MC Cartesian')
py.plot(X, err(QMC_sph), label='QMC Spherical')
py.plot(X, err(QMC_car), label='QMC Cartesian')
py.plot(X, c*X**-.5, 'k--', X, c*X**-1, 'k--', X, c*X**-1.5, 'k--')
py.xlim(1e5, 1e7)
py.ylim(1e-11, 1e0)
py.legend(ncol=2)
py.loglog()
py.show()

