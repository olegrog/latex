#!/usr/bin/env python

import numpy as np
import pylab as py

b1, N, R = 0.1, 16, 8.
p=2

A = (1.+p)/N**(1+p)*(R-N*b1)
n = np.arange(N)
b = b1 + A*n**p
X = np.append(0, np.cumsum(b))
print(A, X)

f = lambda X, T: X**0*np.exp(-X**2/T)/T**.5
py.plot(X, f(X,1), 'r-*', X, f(X,5), 'g-*')
py.show()

f = lambda X, T: X**1*np.exp(-X**2/T)/T**.5
py.plot(X, f(X,1), 'r-*', X, f(X,5), 'g-*')
py.show()

f = lambda X, T: X**2*np.exp(-X**2/T)/T**.5
py.plot(X, f(X,1), 'r-*', X, f(X,5), 'g-*')
py.show()

