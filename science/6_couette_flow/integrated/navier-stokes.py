#!/usr/bin/python

import pylab as py
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad
from functools import partial

L = .5
N = 500
s = 0.5
deltaU = 1e0
tau = 1e-5
Tmax = 1000
eps = 1e-3
err = 1e-18

gamma1 = 1.270042427
gamma2 = 1.922284066
alpha = 1.25*gamma2/gamma1
h = L/N
X = np.linspace(0, L, N)

def der(T):
    derT = np.zeros(N)
    derT[0] = -(3*T[0] - 4*T[1] + T[2])/2/h
    derT[N-1] = (3*T[N-1] - 4*T[N-2] + T[N-3])/2/h
    derT[1:N-1] = 0.5*(T[2:N] - T[0:N-2])/h
    #derT[0] = (T[1] - T[0])/h
    #derT[1:N] = (T[1:N] - T[0:N-1])/h
    return derT

# initial approximation
T = 1 - deltaU**2/2/alpha*(X**2-0.25)
A = -deltaU
sum_q, sum_u = 0., 0.
T_init = T.copy()
T_old = T.copy()

for i in xrange(Tmax):
    if tau < err:
        break
    T1 = der(T)
    T1[0] = 0.
    T1[N-1] = -deltaU*A/2/alpha
    T2 = der(T1)
    #print T, T1, T2
    #py.plot(X,T1)
    #py.show()
    T3 = der(T2)
    F = s*(2*s-1)*T**(2*s-2)*T1**2 + 4*s*T**(2*s-1)*T1*T2 + T**(2*s)*T3
    mask = F < 0
    ones = np.ones(F.size)
    ones[mask] = 1. #-1e-2
    F = F*ones
    #print F
    if np.all(F*tau >= -eps) and np.all(F*tau <= 2.+eps):
        T_old = T.copy()
        f = (A**2/alpha + s*T**(2*s-1)*T1**2 + T**(2*s)*T2)*ones
        # calculate U and Q
        U = -alpha*T1/A
        Q = -alpha*gamma1*np.sqrt(T)*T1
        # update T and A
        dA = np.trapz(f,X)/L
        if dA < alpha:
            A += np.sqrt(1-dA/alpha) - 1
        T -= tau*f
        T[N-1] = 1
        print "err = %.4e, A = %.4e, U = %.4e, Q = %.4e" % (np.trapz(np.abs(f),X), A, np.trapz(U,X), np.trapz(Q,X))
        #print F
    else:
        T = T_old
        tau /= 1.2
        print "Tau =", tau

py.plot(X,T,X,T_init)
#py.plot(X,T-T_init)
py.show()
#np.savetxt(sys.stdout, np.transpose((Kn, U, Q)), fmt='%1.5e') 
