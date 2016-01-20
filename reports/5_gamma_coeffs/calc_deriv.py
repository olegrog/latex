#!/usr/bin/env python

import sys, os
import numpy as np
import pylab as py

#f(x) = C1 + C2*x**2 + C4*x**4 + C6*x**6
nu = np.array([2., 2./3, -1./15, 1./105])
preK1, preK2 = lambda y: 4*y/3, lambda y: 4*y/5
K1 = np.array([ lambda y: 2 + y**2, lambda y: 11./5+0*y, lambda y: 8./7 + 0*y, lambda y: 16./45 + 0*y ])
K2 = np.array([ lambda y: 2 - y**2, lambda y: -16*y**2/21 + 23./7, lambda y: -16*y**2/63 + 40./21, lambda y: -128*y**2/2079 + 64./99 ])
ih_A = -2**1.5*np.array([-2.5, 1., 0, 0])
ih_B = -2**1.5*np.array([2., 0, 0, 0])
corr_A = np.array([ lambda C: 4./5*C[0], lambda C: 29*C[0]/35 + 2*C[1]/7, lambda C: 43*C[0]/189 + 55*C[1]/189 + 4*C[2]/27 ])
corr_B = np.array([ lambda C: 4./7*C[0], lambda C: 5*C[0]/7 + 2*C[1]/9, lambda C: 67*C[0]/297 + 25*C[1]/99 + 4*C[2]/33 ])

def interp(data, deg=16):
    return np.poly1d(np.polyfit(eta, data, deg))

def integrate(func):
    return np.trapz(func(xi), xi)

def calc_coeffs(func, pre, kernel, ih, corr):
    integral = lambda n: integrate(lambda y: pre(y)*kernel[n](y)*func(y)*np.exp(-y**2))
    C = np.zeros(nu.size)
    for i in xrange(nu.size):
        if i:
            ih[i] += corr[i-1](C)
        C[i] = (integral(i) - ih[i] - np.sum(nu[i:0:-1]*C[:i]))/nu[0]
    return C

k = lambda i: i.reshape(-1,1)
asym = lambda c, x, n: np.sum(np.fromfunction(lambda i: c[k(i)]*x**(2*k(i)), (n,), dtype=int), axis=0)
der_asym = lambda c, x, n: np.sum(np.fromfunction(lambda i: (2*k(i)+2)*c[k(i)+1]*x**(2*k(i)+1), (n,), dtype=int), axis=0)

if len(sys.argv) == 1:
    N, xi_cut = int(2e3), 5.4
    xi = np.linspace(0, xi_cut, N+1)
    eta, A_, B_, _, _, _ = np.loadtxt("./tables/ci-functions.txt").T
    A, B = interp(A_), interp(B_)
    derA, derB = A.deriv(), B.deriv()
else:
    dirname = sys.argv[1]
    xi, A_ = np.loadtxt(os.path.join(dirname, 'A.txt')).T
    xi, B_ = np.loadtxt(os.path.join(dirname, 'B.txt')).T
    A, B = lambda xi: A_, lambda xi: B_
    h = xi[1] - xi[0]
    derA_, derB_ = np.gradient(A_, h, edge_order=2), np.gradient(B_, h, edge_order=2)
    derA, derB = lambda xi: derA_, lambda xi: derB_

A_coeffs = calc_coeffs(A, preK1, K1, ih_A, corr_A)
B_coeffs = calc_coeffs(B, preK2, K2, ih_B, corr_B)
np.savetxt(sys.stderr, (A_coeffs, B_coeffs), fmt='%.8f')

def compare_plot(func, coeffs, der):
    print func(xi).shape, asym(coeffs, xi, 2).shape
    py.plot(xi, func(xi)-asym(coeffs, xi, 2), 'r-', xi, func(xi)-asym(coeffs, xi, 3), 'g-', xi, func(xi)-asym(coeffs, xi, 4), 'b-')
    py.axhline(0, c='k', ls=':', lw=.5)
    py.show()
    with np.errstate(divide='ignore'):
        py.plot(xi, der(xi)/xi, 'r-', xi, der_asym(coeffs, xi, 3)/xi, 'g-')
    py.show()

mask = xi<1
xi = xi[mask]
if len(sys.argv) > 1:
    A_, B_ = A_[mask], B_[mask]
    derA_, derB_ = derA_[mask], derB_[mask]

compare_plot(A, A_coeffs, derA)
compare_plot(B, B_coeffs, derB)
