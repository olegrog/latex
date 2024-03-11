#!/usr/bin/env python3

import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, log, exp
from scipy import optimize
sys.path.append('../../09_bgk_lattice/couette')
import vgrid

parser = argparse.ArgumentParser(description='Solver for the Boltzmann--Shakhov equation')
parser.add_argument('-N', type=int, default=20, help='number of points along xi_r')
parser.add_argument('-t', '--timestep', type=float, default=0.2, help='time step')
parser.add_argument('-e', '--end', type=float, default=10, help='total time')
parser.add_argument('--tau', type=float, default=1, help='relaxation time')
parser.add_argument('-f', '--function', default='piecewise', help='type of VDF')
parser.add_argument('-r', '--ratio', type=float, default=2, help='ratio of values in VDF')
parser.add_argument('-c', '--cut', type=float, default=2, help='cutting radius of VDF')
parser.add_argument('-R', '--radius', type=float, default=5, help='radius of the velocity grid')
parser.add_argument('-g', '--grid', default='uniform', metavar='uniform|hermite|polynomial|geometric', help='type of the grid')
parser.add_argument('-p', '--plot', type=int, default=10, help='plot every <int> steps')
parser.add_argument('--gratio', type=float, default=1.15, help='ratio for the geometric grid')
parser.add_argument('--w-min', type=float, default=0.1, help='minimum weight for the polynomial grid')
parser.add_argument('--poly', type=float, default=2, help='power for the polynomial grid')
parser.add_argument('--pdf', action='store_true', help='save PDF file instead')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

### Constants
class fixed:
    D = 1                   # dimension of velocity space
    Pr = 2/3                # Prandtl number

sqr = lambda x: x*x
_c = lambda vel: xi - vel
_cc = lambda vel: sqr(xi - vel)

def calc_macro(f):
    rho = np.sum(f)
    speed = np.sum(f*xi)/rho
    c, cc = _c(speed), _cc(speed)
    temp = 2*np.sum(f*cc)/rho
    qflow = np.sum(f*c*cc)
    hfunc = np.sum(f*log(f/w))
    return rho, speed, temp, qflow, hfunc

def maxwell():
    ones = np.ones_like(w)
    maxw5 = lambda m: w * m[0] * (np.pi*m[2])**(-fixed.D/2) * np.exp(-_cc(m[1]) / m[2])
    maxw5_jac = lambda m: np.vstack((
        ones / m[0],
        2 * _c(m[1]) / m[2],
        (_cc(m[1]) / m[2] - fixed.D/2) / m[2]
    ))

    psi = np.vstack(( ones, xi, sqr(xi) ))
    macro = np.array([ rho, speed, temp ])
    moments = rho*np.array([ 1, speed, fixed.D*temp/2 ])
    vdf = maxw5
    vdf_jac = lambda m: np.einsum('i,li->li', maxw5(m), maxw5_jac(m))

    sol = optimize.root(
        fun = lambda m, psi, moments: np.einsum('li,i', psi, vdf(m)) - moments,
        jac = lambda m, psi, _: np.einsum('li,mi', psi, vdf_jac(m)),
        x0 = macro, args=(psi, moments)
    )
    return vdf(sol.x)

def initial_vdf():
    f = np.zeros_like(xi)
    if args.function == 'piecewise':
        c = args.cut
        if (c > args.radius or c/args.ratio > args.radius):
            raise ValueError('Cutting radius is too big!')
        neg, pos = (-c < xi) * (xi < 0), (0 < xi) * (xi < c/sqrt(args.ratio))
        ratio = -np.sum((xi*w)[neg])/np.sum((xi*w)[pos])
        mass = np.sum(w[neg]) + ratio*np.sum(w[pos])
        f[neg] = 1/mass
        f[pos] = ratio/mass
        if args.verbose:
            print(f'VDF: left = {f[neg][0]:.3g}, right = {f[pos][0]:.3g}, ratio = {ratio:.3g}')
    else:
        raise ValueError('Function type is undefined!')
    return f*w

def print_macro(f):
    rho_, speed_, temp_, qflow_, hfunc_ = calc_macro(f)
    print('Macroscopic variables:')
    print(f' -- density = {rho_:.5g}')
    print(f' -- velocity = {speed_:.5g}')
    print(f' -- temperature = {temp_:.5g}')
    print(f' -- heat flux = {qflow_:.5g}')
    print(f' -- H-function = {hfunc_:.5g}')

####################################################################################################

vgrid_params = {
    'polynomial': { 'w_min': args.w_min, 'p': args.poly },
    'geometric': { 'q': args.ratio },
}
grid = getattr(vgrid, args.grid.capitalize())(args.radius, args.N, **vgrid_params.get(args.grid, {}))
xi, w = grid.x, grid.w
f_0 = initial_vdf()
rho, speed, temp, qflow0, h0 = calc_macro(f_0)
f_M = maxwell()
_, _, _, _, h_M = calc_macro(f_M)
Q1, Q2, H, X = [], [], [], []

for i in range(int(args.end // args.timestep)):
    time = i*args.timestep
    if args.verbose:
        print(f'# Iteration = {i}, time = {time:.5g}')
    e1, e23 = exp(-time/args.tau), exp(-2/fixed.D*time/args.tau)
    q = qflow0*e23
    A = (1-fixed.Pr)/5/pi**(fixed.D/2)*4*q/rho/temp**2*_c(speed)*(_cc(speed)/temp - 5/2)
    f = f_0*e1 + f_M*(1 - e1 + fixed.D*A*(1 - e23))

    Q1.append(q)
    rho, speed, temp, qflow, hfunc = calc_macro(f)
    Q2.append(qflow)
    H.append(hfunc)
    X.append(time)

    if args.verbose:
        print(f'Heat flux = {q}')
        print_macro(f)
    if args.plot > 0 and i % args.plot == 0:
        plt.plot(xi, f/w, '-o', xi, f_M/w, '--')
        plt.axhline(lw=.5, c='k', ls=':')
        plt.axvline(lw=.5, c='k', ls=':')
        #plt.semilogy()
        if args.pdf:
            plt.savefig(f'vdf-{time}.pdf', bbox_inches='tight')
            plt.close()
        else:
            plt.show()

X = np.array(X); H = np.array(H)
Q1 = np.array(Q1); Q2 = np.array(Q2)
plt.plot(X, -Q1, label='Q1')
plt.plot(X, -Q2, label='Q2')
plt.title('Heat flux')
plt.legend()
plt.semilogy()
plt.show()

plt.plot(X, H)
plt.axhline(y=h_M, lw=.5, c='k', ls=':')
plt.semilogx()

if args.pdf:
    plt.savefig('solution.pdf', bbox_inches='tight')
else:
    plt.title('H-function')
    plt.show()
