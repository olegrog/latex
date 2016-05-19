#!/usr/bin/env python
r"Script for evaluating the errors in the various moments calculation of the discrete distribution function"
r"Usage: ./discrete_error.py <U> <N_r> <N_y> <cut_r> <ratio>"
r"Example:  ./discrete_error.py .5 12 22 4.3 1.3"

import sys
import numpy as np
from functools import partial
from collections import namedtuple

_, temp_ratio, N_r, N_y, cut_r, q_r, q_y = sys.argv

Rho, Temp, Speed = 1., 1., np.zeros(3)
Speed[0] = 1e-1
Qflow, Tau = np.array([0.1, 0.2, 0]), np.array([0.05, 0., 0.2])*0
temp_ratio = float(temp_ratio)
heat = temp_ratio*Speed[0]
N_r, N_y = int(N_r), int(N_y)
cut_r, cut_y = float(cut_r), float(cut_r)
q_r, q_y = float(q_r), float(q_y)

zeros, ones = np.zeros_like(Speed), np.ones_like(Speed)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1

### Define different grids

Grid = namedtuple('Grid', ['i2h','i2xi'])

grids = {
    'uniform': Grid(
        i2h =  lambda q, cut, N: cut/N + idx(N)*0,
        i2xi = lambda q, cut, N: idx(N)*cut/N
    ),
    'geometric': Grid(
        i2h =  lambda q, cut, N: h1(q, cut, N) * symm_h(q**np.arange(N)),
        i2xi = lambda q, cut, N: h1(q, cut, N) * (1 if q==1 else symm_xi((q**np.arange(N+1)-1)/(q-1)))
    ),
    'hermite': Grid(
        i2h =  lambda q, cut, N: hermite(cut, N)[0],
        i2xi = lambda q, cut, N: hermite(cut, N)[1]
    ),
    'quadratic': Grid(
        i2h =  lambda q, cut, N: quadratic(q, cut, N)[0],
        i2xi = lambda q, cut, N: quadratic(q, cut, N)[1]
    )
}

def hermite(cut, N):
    H_xi, H_h = np.polynomial.hermite.hermgauss(2*N)
    H_h /= np.exp(-H_xi**2)
    C = 2*cut / np.sum(H_h)
    H_xi *= C
    H_h *= C
    return H_h, H_xi

def quadratic(q, cut, N):
    p=2
    if cut-N*q < 0:
       raise NameError('q = %g is too big' % q) 
    A = (cut-N*q) / np.sum(np.arange(N)**p)
    h = q + A*np.arange(N)**p
    X = np.append(0, np.cumsum(h))
    return symm_h(h), symm_xi(X)

symm_h = lambda x: np.hstack((x[::-1], x))
semi_sum = lambda x: .5*(x[:-1] + x[1:])
symm_xi = lambda x: np.hstack((-semi_sum(x)[::-1], semi_sum(x)))
h1 = lambda q, cut, N: cut/N if q==1 else cut*(q-1)/(q**N-1)
idx = lambda N: np.arange(2*N) - N + .5

### Create grid

gtype_r = gtype_y = 'quadratic'
H_r, Xi_r = grids[gtype_r].i2h(q_r, cut_r, N_r), grids[gtype_r].i2xi(q_r, cut_r, N_r)
H_y, Xi_y = grids[gtype_y].i2h(q_y, cut_y, N_y), grids[gtype_y].i2xi(q_y, cut_y, N_y)
E_r, E_y = np.ones_like(Xi_r), np.ones_like(Xi_y)

dxi = np.einsum('i,j,k', H_r, H_y, H_r)
xi = lambda v=zeros: np.einsum('li,lj,lk->ijkl', (Xi_r-v[0], E_r, E_r), (E_y, Xi_y-v[1], E_y), (E_r, E_r, Xi_r-v[2]))
sqr_xi = lambda v=zeros,s=ones: np.linalg.norm(np.einsum('ijkl,l->ijkl', xi(v), s), axis=3)**2
Maxwell = lambda vel, temp, rho=Rho: rho/(np.pi*temp)**1.5 * np.exp(-sqr_xi(vel)/temp)*dxi

radii = np.array([cut_r, cut_y, cut_r])
ball = sqr_xi(zeros, 1./radii) <= 1

print_grid = False
if print_grid:
    print "--- r:", H_r, Xi_r
    print "--- y:", H_y, Xi_y

sizes = lambda H: (H.min(), H.max(), H.max()/H.min())
print "cut_r = %g, cut_y = %g, q_r = %g, q_y = %g" % (cut_r, cut_y, q_r, q_y)
print "Cell size (r): min = %.4g, max = %.4g, ratio = %.3g" % sizes(H_r)
print "Cell size (y): min = %.4g, max = %.4g, ratio = %.3g" % sizes(H_y)
total = lambda X, cut: np.sum(np.abs(X) <= cut)/2
print "Total cells: %d (%d, %d, %d):" % (np.sum(ball), N_r, N_y, N_r)

def splot(f):
    import pylab as py
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = py.figure()
    ax = p3.Axes3D(fig)
    xv, yv = np.meshgrid(Xi_r, Xi_y, sparse=False, indexing='ij')
    f /= dxi
    f[np.invert(ball[:,:,N_z])] = np.NaN
    #ax.plot_wireframe(xv, yv, np.log(f[:,:,N_z]), linewidth=0.25)
    ax.plot_wireframe(xv, yv, f[:,:,N_z], linewidth=0.25)
    py.savefig('tmp.pdf')
    py.show()

err = lambda theor, real: (real-theor)/theor

def calc_macro(f):
    f = np.einsum('ijk,ijk->ijk', f, ball)
    rho = np.einsum('ijk->', f)
    speed = np.einsum('ijk,ijkl', f, xi())/rho
    c, sqr_c = xi(speed), sqr_xi(speed)
    csqr_c = np.einsum('ijk,ijkl->ijkl', sqr_c, c)
    cc_ij = np.einsum('ijkl,ijkm->ijklm', c, c)
    cc = np.einsum('ijklm,lmn', cc_ij, hodge)
    temp = 2./3*np.einsum('ijk,ijk', f, sqr_c)/rho
    qflow = np.einsum('ijk,ijkl', f, csqr_c)
    tau = 2*np.einsum('ijk,ijkl', f, cc)
    return rho, temp, speed, qflow, tau

def test_1(Temp, Speed):
    P = Rho * Temp
    Xi, Sqr_xi = xi(Speed), sqr_xi(Speed)
    f1 = 1./P/Temp * np.einsum('ijkl,ijkm,n,lmn->ijk', Xi, Xi, Tau, hodge)
    f2 = 4./5/P/Temp * (np.einsum('ijkl,ijk,l->ijk', Xi, Sqr_xi/Temp-2.5, Qflow))
    f = Maxwell(Speed, Temp) * (1 + f1 + f2)
    rho, temp, speed, qflow, tau = calc_macro(f)
    #splot(f)
    print "\n-- Test #1: continual flows - Grad's 13-moment approximation (temp = %g, speed = %g)" % (Temp, Speed[0])
    print "rho =", err(Rho, rho), Rho, rho
    print "temp =", err(Temp, temp), Temp, temp
    print "speed =", err(Speed, speed)
    print "qflow =", err(Qflow, qflow/rho)
    print "tau =", err(Tau, tau/rho)

def test_2(Temp):
    Temp1, Temp2 = Temp, Temp*temp_ratio
    double_temp = np.sqrt(Temp1*Temp2)
    ss_temp = np.sqrt(Temp1) + np.sqrt(Temp2)
    delta_temp = Temp2 - Temp1

    Rho1, Rho2 = 2*Rho*np.sqrt(Temp2)/ss_temp, 2*Rho*np.sqrt(Temp1)/ss_temp
    f = Maxwell(Speed, Temp1, Rho1)
    negative = xi()[...,1] < 0
    f[negative] = Maxwell(-Speed, Temp2, Rho2)[negative]
    rho, temp, speed, qflow, tau = calc_macro(f)

    Rho_ = Rho
    Temp_ = double_temp * (1 + 8./3*np.dot(Speed,Speed)/ss_temp**2)
    Speed_ = Speed * delta_temp / ss_temp**2
    for i in [1,2]:
        Speed_[i] = 1.
        speed[i] += 1.
    Qflow_ = [ 0, -2*Rho*double_temp*delta_temp/ss_temp/np.sqrt(np.pi) * (1 + 2*np.dot(Speed,Speed)/ss_temp**2), 0 ]
    Tau_ = [ 0, 0, 4*Rho*Speed[0]*double_temp/ss_temp/np.sqrt(np.pi) ]

    #splot(f)
    #splot(f*c[0]*c[1])

    print "\n-- Test #2: free molecular flows - sum of 2 half-Maxwellians"
    print "rho =", err(Rho_, rho)
    print "temp =", err(Temp_, temp)
    print "speed =", err(Speed_, speed)
    print "qflow =", err(Qflow_, qflow/rho)
    print "tau =", err(Tau_, tau/rho)

def test_3():
    kn = 1
    gamma_1 = 1.270042427
    eta_, A_, B_, D1_, D2_, F_ = np.loadtxt("../tables/ci-functions.txt").T
    B = np.poly1d(np.polyfit(eta_, B_, 16))
    Xi, Sqr_xi = xi(), sqr_xi()
    phi = 2*Speed[0]* Xi[...,0] * (1-Xi[...,1]*B(np.sqrt(Sqr_xi)))
    f = Maxwell(zeros, Temp) * (1 + phi)
    rho, temp, speed, qflow, tau = calc_macro(f)

    print "\n-- Test #3: linear case - asymptotic solution for hard-sphere molecules"
    print "rho =", err(Rho, rho)
    print "speed =", err(Speed, speed)
    print "tau =", err((0, 0, -2*gamma_1*Speed[0]*kn), tau/rho)

with np.errstate(divide='ignore', invalid='ignore'):
    test_1(Temp*(temp_ratio-1), Speed)
    test_1(Temp, Speed)
    test_1(Temp*temp_ratio, Speed)
    #test_1(Temp*np.sqrt(temp_ratio), Speed)
    Speed[1] = 1e-2
    test_2(Temp)
    #test_2(Temp*temp_ratio)
    #test_2(Temp*np.sqrt(temp_ratio))
