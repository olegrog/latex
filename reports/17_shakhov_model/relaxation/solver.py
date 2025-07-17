#!/usr/bin/env python3

import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, log, exp, fabs
from scipy import optimize
sys.path.append('../../09_bgk_lattice/couette')
import vgrid

parser = argparse.ArgumentParser(description='Solver for the Boltzmann--Shakhov equation')
parser.add_argument('-P', '--Pr', type=float, default=2/3, help='Prandtl number')
parser.add_argument('-N', type=int, default=20, help='number of points along xi_r')
parser.add_argument('-t', '--timestep', type=float, default=0.2, help='time step')
parser.add_argument('-e', '--end', type=float, default=10, help='total time')
parser.add_argument('--tau', type=float, default=1, help='relaxation time')
parser.add_argument('-f', '--function', default='piecewise', help='type of VDF')
parser.add_argument('-a', '--value', type=float, default=1, help='constant value in VDF')
parser.add_argument('-r', '--ratio', type=float, default=2, help='ratio of values in VDF')
parser.add_argument('-i', '--inner', type=float, default=1, help='inner cutting radius of VDF')
parser.add_argument('-o', '--outer', type=float, default=2, help='outer cutting radius of VDF')
parser.add_argument('-R', '--radius', type=float, default=5, help='radius of the velocity grid')
parser.add_argument('-q', '--qflow', type=float, default=1, help='initial heat flux')
parser.add_argument('-g', '--grid', default='uniform', metavar='uniform|hermite|polynomial|geometric', help='type of the grid')
parser.add_argument('-p', '--plot', type=int, default=10, help='plot every <int> steps')
parser.add_argument('-l', '--log', action='store_true', help='use log axes')
parser.add_argument('--gratio', type=float, default=1.15, help='ratio for the geometric grid')
parser.add_argument('--w-min', type=float, default=0.1, help='minimum weight for the polynomial grid')
parser.add_argument('--poly', type=float, default=2, help='power for the polynomial grid')
parser.add_argument('--pdf', action='store_true', help='save PDF file instead')
parser.add_argument('--non-conservative', action='store_true', help='use a non-conservative scheme')
parser.add_argument('-c', '--convolution', default='none', help='type of convolution used for H-function')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
parser.add_argument('-n', '--negative', action='store_true', help='print information on negative values')
args = parser.parse_args()

### Constants
class fixed:
    D = 1                   # dimension of velocity space

sqr = lambda x: x*x
_c = lambda vel: xi - vel
_cc = lambda vel: sqr(xi - vel)

def calc_macro(f):
    rho = np.sum(f)
    speed = np.sum(f*xi)/rho
    c, cc = _c(speed), _cc(speed)
    temp = 2*np.sum(f*cc)/rho/fixed.D
    qflow = np.sum(f*c*cc)
    return rho, speed, temp, qflow

def calc_hfunc(f, f_M, time = -1):
    pos, neg = f>0, f<0
    hfunc = np.sum(f[pos]*log((f/w)[pos]))
    #hfunc = np.sum(f*log(fabs(f/w)))
    #hfunc = np.sum(f*log((f/w+1)))
    if args.convolution == 'none':
        hfunc2 = hfunc
    else:
        match args.convolution:
            case 'fast':
                F = np.convolve(f, f_M/w, mode='same')
            case 'full':
                F = np.einsum('i,ij', f, F_M)
            case 'const':
                F = f+w
            case _:
                raise ValueError(f'Convolution type {args.convolution} is undefined!')
        Pos, Neg = F>0, F<0
        hfunc2 = np.sum(F[Pos]*log((F/w)[Pos]))
    if time > -1 and args.negative:
        if args.verbose:
            if np.count_nonzero(neg) > 0:
                print(f'Negative VDF at {time}: xi = {xi[neg]}, f = {f[neg]}')
        else:
            if np.count_nonzero(neg) > 0:
                idx1, idxN  = xi[neg][0], xi[neg][-1]
                print(f'Negative VDF at {time}: xi = [{idx1}:{idxN}], min(f) = {np.min(f)}')
            if args.convolution != 'none' and np.count_nonzero(Neg) > 0:
                Idx1, IdxN  = xi[Neg][0], xi[Neg][-1]
                print(f'Negative VDF2 at {time}: xi = [{idx1}:{idxN}], min(f) = {np.min(F)}')
    return hfunc, hfunc2

def maxwell(macro):
    ones = np.ones_like(w)
    maxw5 = lambda m: w * m[0] * (pi*m[2])**(-fixed.D/2) * np.exp(-_cc(m[1]) / m[2])
    maxw5_jac = lambda m: np.vstack((
        ones / m[0],
        2 * _c(m[1]) / m[2],
        (_cc(m[1]) / m[2] - fixed.D/2) / m[2]
    ))

    psi = np.vstack(( ones, xi, sqr(xi) ))
    rho_, speed_, temp_ = macro[:3]
    moments = rho_*np.array([ 1, speed_, fixed.D*temp_/2 + speed_**2 ])
    vdf = maxw5
    vdf_jac = lambda m: np.einsum('i,li->li', maxw5(m), maxw5_jac(m))

    if args.non_conservative:
        return vdf(macro)

    sol = optimize.root(
        fun = lambda m, psi, moments: np.einsum('li,i', psi, vdf(m)) - moments,
        jac = lambda m, psi, _: np.einsum('li,mi', psi, vdf_jac(m)),
        x0 = macro[:3], args=(psi, moments)
    )
    if not sol.success:
        raise Exception(sol)

    return vdf(sol.x)

def grad13(macro):
    qflow = macro[3]
    qc = lambda m: 2 * qflow * _c(m[1])
    cc_5T = lambda m: 2 * _cc(m[1]) / (fixed.D+2) / m[2]
    g_qflow = lambda m: qc(m) * (cc_5T(m) - 1)
    vdf = lambda m: maxwell(macro) * (1 + g_qflow(m) / m[0] / m[2]**2)
    return vdf(macro)

pmaxw_delta = lambda qflow: -sqrt(pi)*qflow/sqrt(4 + pi*qflow**2)

def initial_vdf():
    f = np.zeros_like(xi)
    skip = 0
    for function_type in args.function.split(','):
        match function_type:
            case 'constant':
                f += args.value*w
                skip += 1
            case 'maxwell':
                f += maxwell(np.array([1, 0, 1]))
            case 'grad13':
                f += grad13(np.array([1, 0, 1, args.qflow]))
            case 'piecewise':
                i, o = args.inner, args.outer
                ratio, R = args.ratio, args.radius
                if (i >= o):
                    raise ValueError('Inner radius is larger than the outer one!')
                if (o > R or o/sqrt(ratio) > R):
                    raise ValueError('Cutting radius is too big!')
                neg, pos = (-o < xi) * (xi < -i), (i/sqrt(ratio) < xi) * (xi < o/sqrt(ratio))
                ratio = -np.sum((xi*w)[neg])/np.sum((xi*w)[pos])
                mass = np.sum(w[neg]) + ratio*np.sum(w[pos])
                f0 = np.zeros_like(f)
                f0[neg] = 1/mass
                f0[pos] = ratio/mass
                f += f0*w
                if args.verbose:
                    print(f'VDF: left = {f0[neg][0]:.3g}, right = {f0[pos][0]:.3g}, ratio = {ratio:.3g}')
            case 'pmaxw':
                neg, pos = xi < 0, xi > 0
                delta = pmaxw_delta(args.qflow)
                rho1, rho2 = 1 - delta, 1 + delta
                temp1, temp2 = rho2/rho1, rho1/rho2
                if args.verbose:
                    print(f"delta = {delta}, temp1 = {temp1}, temp2 = {temp2}")
                f[neg] += maxwell(np.array([rho1, 0, temp1]))[neg]
                f[pos] += maxwell(np.array([rho2, 0, temp2]))[pos]
            case _:
                raise ValueError(f'Function type {function_type} is undefined!')
    return f/(len(args.function.split(',')) - skip)

def print_macro(f):
    rho_, speed_, temp_, qflow_ = calc_macro(f)
    print('Macroscopic variables:')
    print(f' -- density = {rho_:.5g}')
    print(f' -- velocity = {speed_:.5g}')
    print(f' -- temperature = {temp_:.5g}')
    print(f' -- heat flux = {qflow_:.5g}')

####################################################################################################

vgrid_params = {
    'polynomial': { 'w_min': args.w_min, 'p': args.poly },
    'geometric': { 'q': args.ratio },
}
grid = getattr(vgrid, args.grid.capitalize())(args.radius, args.N, **vgrid_params.get(args.grid, {}))
xi, w = grid.x, grid.w
f_0 = initial_vdf()
rho, speed, temp, qflow0 = calc_macro(f_0)
macro = np.array([ rho, speed, temp ])
if args.convolution == 'full':
    Xi1 = np.einsum('i,j->ij', xi, np.ones_like(xi))
    Xi2 = np.einsum('i,j->ji', xi, np.ones_like(xi))
    F_M = rho*(pi*temp)**(-fixed.D/2)*np.exp(-(Xi1-Xi2-speed)**2/temp)
f_M = maxwell(macro)
h_0, h2_0 = calc_hfunc(f_0, f_M)
macro = np.array([ rho, speed, temp ])
h_M, h2_M = calc_hfunc(f_M, f_M)
X, Q1, Q2, H, H2 = [], [], [], [], []

### 1. Plot VDF
for i in range(int(args.end // args.timestep)):
    time = i*args.timestep
    if args.verbose:
        print(f'# Iteration = {i}, time = {time:.5g}')
    exp1, expPr = exp(-time/args.tau), exp(-args.Pr*time/args.tau)
    A = 2*qflow0/rho/temp**2*_c(speed)*(2*_cc(speed)/(fixed.D+2)/temp - 1)
    f = f_0*exp1 + f_M*(1 - exp1 + A*(expPr - exp1))
    q = qflow0*expPr

    Q1.append(q)
    rho, speed, temp, qflow = calc_macro(f)
    hfunc, hfunc2 = calc_hfunc(f, f_M, time)
    Q2.append(qflow)
    H.append(hfunc)
    H2.append(hfunc2)
    X.append(time)

    if args.verbose:
        print(f'Heat flux = {q}')
        if i > 0:
            print(f'dH/dt = {(H[i] - H[i-1])/(X[i] - X[i-1])}')
        print_macro(f)
    if args.plot > 0 and i % args.plot == 0:
        plt.plot(xi, f/w, '-', xi, f_M/w, '--')
        plt.axhline(lw=.5, c='k', ls=':')
        plt.axvline(lw=.5, c='k', ls=':')
        #plt.semilogy()
        if args.pdf:
            plt.savefig(f'vdf-{time}.pdf', bbox_inches='tight')
            plt.close()
        else:
            plt.show()

### 2. Plot heat flux
X = np.array(X); H = np.array(H); H2 = np.array(H2)
Q1 = np.array(Q1); Q2 = np.array(Q2)
sgn = np.sign(Q1[0])
plt.plot(X, sgn*Q1, label='exact')
plt.plot(X, sgn*Q2, '--', label='computed')
plt.title('Heat flux')
plt.legend()
plt.semilogy()
if args.pdf:
    plt.savefig(f'heat_flux.pdf', bbox_inches='tight')
    plt.close()
else:
    plt.show()

### 3. Plot H-function
plt.margins(x=0)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.axhline(y=h_M, ls='--', color=colors[1])
plt.plot(X, H)
if args.convolution != 'none':
    plt.plot(X, (H2-h2_M)*(h_0-h_M)/(h2_0-h2_M)+h_M, color=colors[2])
if args.verbose and args.function == 'pmaxw':
    d = pmaxw_delta(args.qflow)
    dd, r = sqr(d), (1 + d)/(1 - d)
    H_0 = log((1-dd)/pi)/2 + d*log(r) - 1/2
    dHdt_BGK = -d*log(r) - dd/(1-dd)
    dHdt_Shakhov = 4/3*(1-args.Pr)/pi*d/sqrt(1-dd)*(log(r) + 2*d/(1-dd))
    print(f"H_0 = {H_0}, dH/dt = {dHdt_BGK + dHdt_Shakhov}")
    tangent = lambda x: H_0 + (dHdt_BGK + dHdt_Shakhov)*x
    tangent2 = lambda x: H[0] + (H[1] - H[0])/(X[1] - X[0])*x
    plt.plot(X, tangent(X), lw=.5, ls='--', c='k')
    plt.plot(X, tangent2(X), lw=.5, ls='--', c='b')
if args.log:
    plt.semilogx()
else:
    plt.axvline(lw=.5, c='k', ls='-')

print('H-function always decreases:', np.all(H[1:] < H[:-1]))
if args.pdf:
    plt.savefig('solution.pdf', bbox_inches='tight')
else:
    plt.title('H-function')
    plt.show()
