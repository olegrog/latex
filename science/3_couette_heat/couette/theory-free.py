#!/usr/bin/python

import pylab as py
import numpy as np
import sys
from scipy.integrate import quad
from functools import partial

a = 2*np.sqrt(np.pi)
f = lambda n,x,t: np.exp(-t*t - x/t)*t**n
T0 = lambda x: quad(partial(f, 0, x), 0, np.infty)[0]
T2 = lambda x: quad(partial(f, 2, x), 0, np.infty)[0]

L = .5

Kn = np.logspace(np.log10(5e-1), np.log10(1e3), num=50)
K = Kn*np.sqrt(np.pi)/2

def quadvec(T):
    f = lambda k,x: T((L-x)/k) - T((L+x)/k)
    return np.vectorize(lambda k: quad(partial(f,k), 0, L)[0]/a)

U = quadvec(T0)(K)
Q = quadvec(T2)(K) - U/2

np.savetxt(sys.stdout, np.transpose((Kn, U, Q)), fmt='%1.4e') 
