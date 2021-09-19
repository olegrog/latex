#!/usr/bin/python

import numpy as np
import pylab as py
import sys
from functools import partial
from scipy.integrate import quad, dblquad
from scipy.interpolate import interp1d
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve

X = np.logspace(-3, 1.8, 50)                           # precision 1e-6
a, b = np.sqrt(np.pi), 1.5*np.euler_gamma
c = np.exp(b)
j = lambda n,x,t: np.exp(-t*t - x/t)*t**n               # Abramowitz function \int j dt
g0 = np.vectorize(lambda x: a/2 + (0 if x==0 else x/(1+x) * (np.log(x) + b-1)))    # regularize T_{0}
g_1 = np.vectorize(lambda x: -(0 if x==0 else np.log(x/(1+np.exp(b)*x))) - b)       # regularize T_{-1}

def interp(f, F0):
    quadvec = lambda ff: np.vectorize(lambda x: np.log(quad(partial(ff,x), 0, np.infty)[0]))
    func = interp1d(np.append(0, X), np.append(np.log(F0), quadvec(f)(X)), fill_value=-np.infty, bounds_error=False)
    return lambda x: np.exp(func(x))

T2 = interp(partial(j,2), a/4)
T1 = interp(partial(j,1), .5)
T0r = interp(lambda x,t: j(0,x,t) / g0(x), 1)
T0 = lambda x: T0r(x) * g0(x)
T_1r = interp(lambda x,t: j(-1,x,t) / g_1(x), 1)
T_1 = lambda x: T_1r(x) * g_1(x)

def plot_abra():
    X = np.logspace(-3,3,500)
    py.plot(X, T0(X), X, T1(X), X, T2(X))
    py.loglog()
    py.show()

L = .5
N = 2000
X = np.linspace(0, L, N)

def printt(Y):
    print("--------------------------")
    np.savetxt(sys.stdout, Y, fmt='%+.2e')
    print("--------------------------")

def splot(X, K):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = py.figure()
    ax = p3.Axes3D(fig)
    xv, yv = np.meshgrid(X, X, sparse=False, indexing='ij')
    ax.plot_wireframe(xv, yv, K, rstride=N/10, cstride=N/10)

def weight_matrix(alpha, beta):
    N = len(X)
    C = np.eye(N)
    for i in range(N):
        C[i, 1:N-1] = alpha(i, np.arange(2,N)) + beta(i, np.arange(1,N-1))
        C[i, 0] = alpha(i, 1)
        C[i, N-1] = beta(i, N-1)
    return C

def apply_abs(i,j,F,G):
    h = L/len(X)
    if (i+j >= 1):
        return F(h*(i+j-1)) - F(h*(i+j))
    elif (i+j <= 0):
        return G(-h*(i+j)) - G(-h*(i+j-1))
    else:
        return 2*F(0) - F(h*(i+j)) - G(-h*(i+j-1))

def solve_linalg(k, T, F0, F1, f):
    N, h = len(X), L/len(X)
    I = np.eye(N)
    S,Y = np.meshgrid(X,X)
    abs_func = np.vectorize(apply_abs)
    F0, F1 = partial(F0, k), partial(F1, k)
    G0 = lambda i,j: abs_func(i, j, F0, F0) - b*h
    G1 = lambda i,j: abs_func(i, j, F1, lambda x: -F1(x)) - b*h*h*(i+j-.5)
    A = weight_matrix(
        lambda i,j: (j-i)*G0(-i,j) - G1(-i,j)/h,
        lambda i,j: G1(-i,j)/h - (j-i-1)*G0(-i,j)
    )
    B = weight_matrix(
        lambda i,j: (j+i)*G0(i,j) - G1(i,j)/h,
        lambda i,j: G1(i,j)/h - (j+i-1)*G0(i,j)
    )
    #splot(X, A*T(np.abs(S-Y)/k)/k)
    #splot(X, B*T((S+Y)/k)/k)
    #py.show()
    phi = solve(a*I - A*T(np.abs(S-Y)/k)/k + B*T((S+Y)/k)/k, f(X))
    p_xy = -(k*(T2(0)-T2(1./k)) + np.trapz((T1((L-X)/k) - T1((L+X)/k))*phi, X))*2/a
    Phi = np.outer(phi,np.ones(N))
    Q = np.trapz(T2((L-X)/k)-T2((L+X)/k) + np.trapz((T1(np.abs(S-Y)/k) - T1((S+Y)/k))*Phi, X)/k, X)/2/a
    #splot(X, K(*XX))
    #py.plot(X, phi)
    w = np.vectorize(lambda x: 0 if x==0 else x*np.log(x)/a)
    ww = lambda x: k*x*x*(2*np.log(x)-1)/4/a
    #print >> sys.stderr, k, np.trapz(phi, X), np.trapz(phi - w((L-X)/k), X) + ww(L/k)
    print(k, p_xy, np.trapz(phi, X)/2, Q, file=sys.stderr)
    #np.savetxt(sys.stdout, np.transpose((X, phi)), fmt='%1.4e')
    return k, p_xy, np.trapz(phi, X)/2, Q

Kn = np.logspace(np.log10(1e-2), np.log10(1e+2), num=40)
Kn = Kn[14:16]
#Kn = np.array([ 1e-1 ])
K = Kn*np.sqrt(np.pi)/2

f = lambda x: T0((L-x)/k) - T0((L+x)/k)
r0 = lambda x: (0 if x==0 else x*(np.log(x) - 1))
r1 = lambda x: (0 if x==0 else .5*x*x*(np.log(x) - .5))
F0 = lambda k,x: r0(x) - r0(k+c*x)/c
F1 = lambda k,x: r1(x) - (r1(k+c*x) - k*r0(k+c*x))/c/c
data = np.array([ solve_linalg(k, T_1r, F0, F1, f) for k in K ])

np.savetxt(sys.stdout, data, fmt='%1.5e') 
