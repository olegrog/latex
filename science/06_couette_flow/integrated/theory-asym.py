#!/usr/bin/python

import pylab as py
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad
from functools import partial

L = .5

xlogx = np.vectorize(lambda x: 0 if x==0 else x*np.log(x))
Fxlogx = lambda x: x*x*(2*np.log(x)-1)/4

def get_consts(kind):
    if kind == "hs":
        g1, g2 = 1.270042427, 1.922284066
        k_0, d_1, = -1.2540, 2.4001
        eta, Y_0, Y_1, Omega_1, Theta_1, H_A, H_B, Y_a4 = py.loadtxt('../tables/kn-layer-hs.txt').T
        Y_1 = 2*Y_1
        return g1, g2, k_0, d_1, eta, Y_0, H_A, Theta_1, 0.43, 0.078, 0.7
    elif kind == "bkw":
        g1, g2 = 1., 1.
        k_0, d_1 = -1.01619, 1.30272
        eta, Y_0 = py.loadtxt('../tables/kn-layer-bkw.txt').T
        H_A = Y_0 / 2
        return g1, g2, k_0, d_1, eta, Y_0, H_A, Theta_1, 0.345, 0.172
    else:
        raise "Invalid potential!"

def print_g(Kn, N):
    X = np.linspace(0, L, N)
    np.savetxt(sys.stdout, np.transpose((X, (.5-2*v_x(Kn*np.sqrt(np.pi)/2,X))[::-1])), fmt='%1.4e')

def interp(eta, F):
    log_interp = interp1d(eta, np.log(F), kind='linear')
    return lambda x: np.exp(log_interp(x))

def quad_log(F, s, x):
    a = min(1,x)
    b = min(25,x)
    return quad(lambda t: F(t)-s*xlogx(t), 0, a)[0] + s*Fxlogx(a) + quad(F,a,b)[0]

g1, g2, k_0, d_1, eta, Y_0, H_A, Theta_1, sY, sH, sTheta = get_consts(sys.argv[1])
alpha = 1.25*g2/g1
Y_0 = interp(eta, Y_0)
H_A = interp(eta, H_A)
Theta_1 = interp(eta, Theta_1)
Kn = np.logspace(np.log10(2e-3), np.log10(2e0), num=30)
K = Kn*np.sqrt(np.pi)/2

kn_layer = lambda F, s, k: k*k*(2*quad_log(F,s,L/k) - quad_log(F,s,2*L/k))
kn_layer2 = lambda F, s, k: k*k*quad_log(F,s,2*L/k)

#U = [ (.125 + kn_layer(Y_0, sY, k))/(1-2*k_0*k) for k in K ]
U = [ (.125 + kn_layer(Y_0, sY, k))/(1-2*k_0*k) for k in K ]
Q = [ -kn_layer(H_A, sH, k)/(1-2*k_0*k) for k in K ]
tau = [ ( 1./6 + d_1*k - kn_layer2(Theta_1, sTheta, k) )/4/alpha/(1-2*k_0*k)**2 for k in K ]
#tau = [ ( (1+4*k_0*k)/6 + d_1*k - kn_layer2(Theta_1, sTheta, k) )/4/alpha for k in K ]

np.savetxt(sys.stdout, np.transpose((Kn, U, Q, tau)), fmt='%1.5e') 
