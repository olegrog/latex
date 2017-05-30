#!/usr/bin/env python
r"Solver for the plane Couette-flow problem"
r"Example: ./solver.py --solver=bgk -N=10 -M=10"

import argparse
import numpy as np
from functools import partial
from collections import namedtuple

parser = argparse.ArgumentParser(description='2D contour plots of scalar and vector from VTK data')
parser.add_argument('-N', type=int, default=8, metavar='value', help='number of cells in the physical space')
parser.add_argument('-M', type=int, default=8, metavar='value', help='number of points along each axis in the velocity space')
parser.add_argument('--time', type=int, default=int(1e2), metavar='value', help='Maximum total time')
parser.add_argument('--radius', type=float, default=4., metavar='value', help='radius of the velocity grid')
parser.add_argument('--solver', default='bgk', metavar='value', help='type of solver')
parser.add_argument('--kn', type=float, default=1e-1, metavar='value', help='Knudsen number')
parser.add_argument('--velocity', type=float, default=1e-1, metavar='value', help='velocity of plate')
args = parser.parse_args()

### Constants
zeros, ones, dom = np.zeros(3), np.ones(3), np.ones(args.N)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
X = (np.arange(args.N) + .5) / args.N / 2

### Velocity grid
idx = lambda N: np.arange(2*N) - N + .5
_H, _Xi = args.radius / args.M + idx(args.M)*0, args.radius / args.M * idx(args.M)
H, Xi = np.einsum('i,a', _H, dom), np.einsum('i,a', _Xi, dom)
_I, I, Zeros = np.ones_like(_Xi), np.ones_like(Xi), np.einsum('l,a', zeros, dom)

dxi = np.einsum('i,j,k', _H, _H, _H)
_xi = lambda v=zeros: np.einsum('li,lj,lk->ijkl', (_Xi-v[0], _I, _I), (_I, _Xi-v[1], _I), (_I, _I, _Xi-v[2]))
_sqr_xi = lambda vel, s=ones: np.linalg.norm(np.einsum('ijkl,l->ijkl', _xi(vel), s), axis=3)**2
_vel = lambda v, i: np.einsum('a,i', v[:,i], _I) 
_Maxwell = lambda rho, vel, temp: rho/(np.pi*temp)**1.5 * np.exp(-_sqr_xi(vel)/temp)*dxi

xi = lambda v=Zeros: np.einsum('lai,laj,lak->aijkl', (Xi-_vel(v,0), I, I), (I, Xi-_vel(v,1), I), (I, I, Xi-_vel(v,2)))
sqr_xi = lambda vel, s=ones: np.linalg.norm(np.einsum('aijkl,l->aijkl', xi(vel), s), axis=4)**2
Maxwell = lambda rho, vel, temp: np.einsum('a,aijk,ijk->aijk', rho/(np.pi*temp)**1.5, np.exp(np.einsum('aijk,a->aijk', -sqr_xi(vel), 1./temp)), dxi)

_gr1 = lambda vel, temp, tau: 2 * np.einsum('ijkl,ijkm,n,lmn->ijk', xi(vel), xi(vel), tau, hodge)
_gr2 = lambda vel, temp, qflow: 4./5 * (np.einsum('ijkl,ijk,l->ijk', xi(vel), sqr_xi(vel)/temp-2.5, qflow))
Grad13 = lambda rho, vel, temp, tau, qflow: Maxwell(vel, temp, rho) * (1 + (_gr1(vel, temp, tau) + _gr2(vel, temp, qflow))/(rho*temp**2))

# Masks
_ball = _sqr_xi(zeros, ones/args.radius) <= 1

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

def plot_profiles(f):
    rho, temp, vel, qflow, tau = calc_macro(f)
    import pylab as py
    py.plot(X, vel[:,0]/args.velocity, 'r', X, 1e3*(rho-1), 'b*--')
    py.show()

def calc_macro(f):
    #f = np.einsum('aijk,ijk->ijk', f, ball)
    rho = np.einsum('aijk->a', f)
    vel = np.einsum('aijk,ijkl,a->al', f, _xi(), 1./rho)
    c, sqr_c = xi(vel), sqr_xi(vel)
    csqr_c = np.einsum('aijk,aijkl->aijkl', sqr_c, c)
    cc_ij = np.einsum('aijkl,aijkm->aijklm', c, c)
    cc = np.einsum('aijklm,lmn', cc_ij, hodge)
    temp = 2./3*np.einsum('aijk,aijk,a->a', f, sqr_c, 1./rho)
    qflow = np.einsum('aijk,aijkl->al', f, csqr_c)
    tau = 2*np.einsum('aijk,aijkl->al', f, cc)
    return rho, temp, vel, qflow, tau

def total_values(f):
    rho, temp, vel, qflow, tau = calc_macro(f)
    print "Total mass =", sum(rho)/args.N

def test():
    sca_ones, vec_zeros = np.array([1.,]), np.array([zeros,])
    rho, temp, vel, qflow, tau = calc_macro(Maxwell(sca_ones, vec_zeros, sca_ones))
    print rho[0], temp[0], vel[0], qflow[0], tau[0]
    # print calc_macro(Grad13([0,0,1e-1],zeros,[1e-1,0,0],1,1))

#test()

def collisions(f, delta_t):
    rho, temp, vel, qflow, tau = calc_macro(f)
    return f + 2*delta_t/(np.sqrt(np.pi)*args.kn)*np.einsum('a, aijk', rho, Maxwell(rho, vel, temp) - f)

T_B = 1
half_M = -np.sum(_xi()[...,1] * _Maxwell(1., np.array([args.velocity, 0, 0]), T_B) * (_xi()[...,1] < 0))
print half_M

def transfer(f, delta_t):
    new_f, N, T_B = np.copy(f), args.N, 1.
    pos, neg = xi()[1:N][...,1] > 0, xi()[0:N-1][...,1] < 0
    antisymm, diffuse = xi()[0][...,1] > 0, xi()[N-1][...,1] < 0
    # Upwind scheme
    new_f[1:N][pos] -= delta_t*np.einsum('ijk, aijk->aijk', _xi()[...,1], f[1:N] - f[0:N-1])[pos]
    new_f[0:N-1][neg] -= delta_t*np.einsum('ijk, aijk->aijk', _xi()[...,1], f[1:N] - f[0:N-1])[neg]
    # Boundary conditions
    new_f[0][antisymm] = f[0,::-1,::-1][antisymm]
    sigma = np.sum(_xi()[...,1] * f[N-1] * (_xi()[...,1] > 0)) / half_M
    new_f[N-1][diffuse] = _Maxwell(sigma, np.array([args.velocity, 0, 0]), T_B)[diffuse]
    return new_f

def solve_bgk():
    h = .5/args.N
    delta_t = h/args.radius # from CFL
    Rho0, Vel0, Temp0 = dom, np.einsum('a,l', 2*X*args.velocity, [1.,0,0]), dom
    f = Maxwell(Rho0, Vel0, Temp0)
    rho0, temp0, vel0, qflow0, tau0 = calc_macro(f)
    total_values(f)
    for i in xrange(args.time):
        f = transfer(f, delta_t/2/h)
        f = collisions(f, delta_t)
        f = transfer(f, delta_t/2/h)
        if not i % 10:
            rho, temp, vel, qflow, tau = calc_macro(f)
            print "Rho error =", sum((rho-rho0)**2)/args.N
    total_values(f)
    rho, temp, vel, qflow, tau = calc_macro(f)
    print rho, temp, vel, qflow, tau
    plot_profiles(f)
    #splot(f)

with np.errstate(divide='ignore', invalid='ignore'):
    {
        'bgk': solve_bgk,
        'lbe': solve_bgk,
        'hybrid': solve_bgk,
    }[args.solver]()
