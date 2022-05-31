#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.special import erf, erfc
from numpy import sqrt, exp

parser = argparse.ArgumentParser(description='Solve a Stefan problem')

str2pair = lambda s: [float(item) for item in s.split(':')]
str2list = lambda s: [float(item) for item in s.split(',')]

parser.add_argument('-N', type=int, default=100, help='number of points along each axis')
parser.add_argument('-Tm', type=float, default=0, help='melting temperature')
parser.add_argument('-T0', type=float, default=2, help='initial temperature')
parser.add_argument('-Tb', type=float, default=-4, help='boundary temperature')
parser.add_argument('-LbyCs', type=float, default=4, help='latent heat/heat capacity in solid')
parser.add_argument('-K', '--Kratio', type=float, default=1, help='K_S/K_L')
parser.add_argument('-k', '--kratio', type=float, default='1', help='kappa_S/kappa_L')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
parser.add_argument('-d', '--debug', action='store_true', help='maximum information')
parser.add_argument('-o', '--output', type=str, default='profiles.pdf', help='PDF filename')
parser.add_argument('-s', '--figsize', type=str2pair, default='5:4', help='figure size')
parser.add_argument('--xrange', type=str2pair, default='0:10', help='range of values along x-axis')
parser.add_argument('--trange', type=str2pair, default='1e-2:20', help='range of values along t-axis')
parser.add_argument('--xpoints', type=str2list, default='0.5,1,2', help='list of x values')
parser.add_argument('--tpoints', type=str2list, default='0.1,1,10,100', help='list of t values')
parser.add_argument('--pdf', action='store_true', help='save a PDF file instead')
args = parser.parse_args()

class Style:
    thin = { 'linestyle': '-', 'linewidth': 0.5, 'color': 'k' }
    dashed = { 'linestyle': '--', 'linewidth': 0.5, 'color': 'k' }
    point = { 'color': 'black', 'marker': '.', 'linestyle': 'None', 'zorder': 10 }

if args.pdf:
    matplotlib.use('pgf')
    params = {
        'axes.labelsize': 11,
        'font.size': 11,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'figure.figsize': args.figsize
    }
    plt.rcParams.update(params)

_xmesh = lambda x0: np.r_[np.linspace(args.xrange[0], x0, args.N),
    np.geomspace(x0, args.xrange[1], args.N)]
_tmesh = lambda t0: np.fromiter([ t for t in np.r_[np.geomspace(args.trange[0], t0, args.N),
    np.geomspace(t0, args.trange[1], args.N)] if t <= args.trange[1] ], dtype=type(t0))
_front = lambda t: 2*lambda_*sqrt(t)
_front_inv = lambda x: (x/2/lambda_)**2
_lambda_eq = lambda l: exp(-l**2)/erf(l) - sqrt(args.kratio)/args.Kratio*(args.T0 - args.Tm) \
    /(args.Tm - args.Tb)*exp(-l**2*args.kratio)/erfc(l*sqrt(args.kratio)) \
    - l*args.LbyCs*(args.Tm - args.Tb)
lambda_ = root_scalar(_lambda_eq, bracket=[1e-2, 1e1]).root

if args.debug:
    X = np.logspace(-2, 1, args.N)
    plt.plot(X, _lambda_eq(X))
    plt.plot(lambda_, 0, **Style.point)
    plt.semilogx()
    plt.axhline(**Style.dashed)
    plt.show()

if args.verbose:
    print(f'lambda = {lambda_}')

def exact(x, t):
    X_ = _front(t)
    return np.where(x < X_,
        (args.Tm - args.Tb)/erf(lambda_)*erf(lambda_*x/X_) + args.Tb,
        args.T0 - (args.T0 - args.Tm)/erfc(lambda_*sqrt(args.kratio))*erfc(lambda_*x/X_))

fig, axs = plt.subplots(ncols=2, figsize=args.figsize[0]*np.array((2.2,1)))

### 1. Draw the exact solution
for x in args.xpoints:
    T = _tmesh(_front_inv(x))
    axs[0].plot(T, exact(x,T), label=f'$x={x}$')
    if args.verbose and _front_inv(x) < args.trange[1]:
        axs[0].plot(_front_inv(x), args.Tm, **Style.point)

axs[0].axhline(y=args.Tm, **Style.dashed)
axs[0].set_xlabel('$t$')
axs[0].set_ylabel('$T$', rotation=0)
axs[0].legend()
if args.verbose:
    axs[0].plot(0, args.T0, **Style.point)

for t in args.tpoints:
    X = _xmesh(_front(t))
    axs[1].plot(X, exact(X,t), label=f'$t={t}$')
    if args.verbose:
        axs[1].plot(_front(t), args.Tm, **Style.point)

axs[1].axhline(y=args.T0, **Style.dashed)
axs[1].axhline(y=args.Tm, **Style.dashed)
axs[1].set_xlabel('$x$')
axs[1].set_ylabel('$T$', rotation=0)
axs[1].legend()
if args.verbose:
    axs[1].plot(0, args.Tb, **Style.point)

plt.tight_layout(pad=1)
if args.pdf:
    if args.verbose:
        print(f' -- Save to {args.output}')
    plt.savefig(args.output, bbox_inches='tight')
else:
    plt.show()

