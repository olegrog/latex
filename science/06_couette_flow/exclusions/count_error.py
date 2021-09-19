#!/usr/bin/env python

import numpy as np
import pylab as py
import random
from scipy.optimize import curve_fit

X0 = np.array([2])**(-.5)
X0 = np.array([2,3,5,7,11,13,17,19,23])**(-.5)

s = lambda x: (1 + np.sign(x-x0))/2

data = [
    ( 1, lambda x, x0: np.sign(x-x0), lambda x0: 1 - 2*x0 ),
    ( 2, lambda x, x0: np.absolute(x-x0), lambda x0: x0**2 - x0 + .5 ),
    ( 2, lambda x, x0: (x*(x0-1+s(x)) - s(x)*x0)*(x-x0), lambda x0: x0*(-2*x0**2 + 3*x0 - 1)/6 )
#    ( 2, lambda x, x0: x**3, lambda x0: .25 )
]

def exclude(arr, M):
    return np.array(random.sample(arr, np.size(arr)-M))

def integr(f, N, x0, M=0):
    Y = np.linspace(0, 1, N+1)
    XX = .5*(Y[0:N] + Y[1:N+1])
    if M > 0:
        # M = int(N/2)  # --> it gives an error O[N**(-.5)]
        XX = exclude(XX, M)
    return np.sum(f(XX, x0))/(N-M)
    #X = np.linspace(0, 1, N)
    #return np.sum(f(X, x0))/N
    #return np.trapz(f(X, x0), X)

for i, f, J in data:
    for x0 in X0:
        eps = 1e-17
        N, N0 = np.round(np.logspace(3, 5, 50)), 1e5
        M, M0 = np.int64(np.round(np.logspace(1, 3, 20))), int(0.7e3)
        integr_ = np.vectorize(integr)
        func = lambda x, a, b: a*x + b
        p, m = curve_fit(func, np.log(N), np.log(np.abs(integr_(f, N, x0, M0) - J(x0))+eps), [-1, 0])
        q, m = curve_fit(func, np.log(M), np.log(np.abs(integr_(f, N0, x0, M) - J(x0))+eps), [-1, 0])
        print("%.5f %+.5f %.3f %.3f" % (x0, J(x0), p[0], q[0]))

        py.loglog()
        #py.plot(M, np.abs(integr_(f, N0, x0, M) - J(x0)))
        #py.show()
    print()
