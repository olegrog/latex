#!/usr/bin/env python

import sys
import numpy as np
from scipy.integrate import quad, dblquad

def interp(data, deg=16):
    return np.poly1d(np.polyfit(eta, data, deg))

cut = 5.5
eta, A_, B_, _, _, _ = np.loadtxt("./tables/ci-functions.txt").T
A, B = interp(A_), interp(B_)

print 3.82368e-01

zzzz = lambda x: x**3 * np.exp(-x**2) * B(x)
print quad(zzzz, 0, cut)[0]**2 * (2**1.5/15)
