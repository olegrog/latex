#!/usr/bin/python

import numpy as np
import sys
import os
from scipy.interpolate import interp1d
from scipy.integrate import quad
from functools import partial

L = .5

xlogx = np.vectorize(lambda x: 0 if x==0 else x*np.log(x))
Fxlogx = lambda x: x*x*(2*np.log(x)-1)/4

g1, g2 = 1.270042427, 1.922284066
k_0, d_1 = -1.2540, 2.4001
eta, Y_0, Y_1, Omega_1, Theta_1, H_A, H_B, Y_a4 = np.loadtxt('../tables/kn-layer-hs.txt').T
Y_1 = 2*Y_1
sY, sH, sTheta, sOmega = 0.43, 0.077, 0.68, 0.50

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

alpha = 1.25*g2/g1
Y_0 = interp(eta, Y_0)
H_A = interp(eta, H_A)
Theta_1 = interp(eta, Theta_1)
Omega_1 = interp(eta, Omega_1)

kn_layer = np.vectorize(lambda F, s, k: k*k*(2*quad_log(F,s,L/k) - quad_log(F,s,2*L/k)))
kn_layer2 = np.vectorize(lambda F, s, k: k*k*quad_log(F,s,2*L/k))

def expansion(row):
    U, Pxy1, Pxy2, M0, M1, Qx2, Qy1, Qy2, Pxx2, Pyy2, Pzz2, T0, T1, DU0, DT0 = row
    Pxy = K*Pxy1 + K*K*Pxy2
    #Pxy = K*Pxy1 / (1 - K*Pxy2/Pxy1)
    M, T = M0 + K*M1, T0 + K*T1
    Qx, Qy = K**2*Qx2, K*Qy1 + K**2*Qy2
    Pxx, Pyy, Pzz = K**2*np.array([Pxx2, Pyy2, Pzz2])
    return U, Pxy, M, Qx, Qy, Pxx, Pyy, Pzz, T

def knudsen_layer(Kn, DU0, DT0, M, T, Qx, Pxx, Pzz, T_B=1):
    K = Kn*np.sqrt(np.pi)/2
    M += DU0*kn_layer(Y_0, sY, K*T_B)
    T += DT0*kn_layer2(Theta_1, sTheta, K*T_B)
    Qx -= DU0*kn_layer(H_A, sH, K*T_B)
    ThetaPlusOmega = lambda x: Theta_1(x) - Omega_1(x)
    Pxx += 1.5*DT0*kn_layer2(ThetaPlusOmega, sTheta-sOmega, K*T_B)/T_B 
    Pzz += 1.5*DT0*kn_layer2(ThetaPlusOmega, sTheta-sOmega, K*T_B)/T_B

def save_data(filename):
    np.savetxt(filename, np.transpose((Kn, Pxy, M, Qx, Qy, Pxx, Pyy, Pzz, T)),
        fmt='%1.5e', header='       Kn          Pxy           M           Qx          Qy         Pxx         Pyy         Pzz           T')


if sys.argv[1] == 'ns':
    for row in np.loadtxt('./navier-stokes.txt'):
        U, Pxy, M, Qx, Qy, Pxx, Pyy, Pzz, T = expansion(row)
        print U
        Kn = np.logspace(np.log10(1e-3), np.log10(1e0), num=31)
        knudsen_layer(Kn, DU0/U, DT0/U, M, T, Qx, Pxx, Pzz)
        save_data("asym-" + str(U) + ".txt")
elif sys.argv[1] == 'asym':
    for filename in os.listdir('.'):
        if filename.startswith("ns-"):
            print filename
            U = float(filename.split('-')[1].rsplit('.', 1)[0])
            Kn, Pxy, M, Qx, Qy, Pxx, Pyy, Pzz, T, DU0, DT0, T_B = np.loadtxt(filename).T
            knudsen_layer(Kn, DU0/U, DT0/U, M, T, Qx, Pxx, Pzz, T_B)
            save_data(filename.replace('ns', 'asym'))
else:
    raise Exception('Illegal type!')


