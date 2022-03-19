#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root
from scipy.special import erf
from numpy import sqrt, exp, log
from functools import partial

parser = argparse.ArgumentParser(
    description='Draw a solid--liquid interface of a melt pool')

str2pair = lambda s: [float(item) for item in s.split(':')]

parser.add_argument('-N', type=int, default=100, help='number of points along each axis')
parser.add_argument('-T', '--Tratio', type=float, default=0.05, help='T_M/T_0')
parser.add_argument('-s', '--figsize', type=str2pair, default='5:4', help='figure size')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
parser.add_argument('-d', '--debug', action='store_true', help='maximum information')
parser.add_argument('-o', '--output', type=str, default='profiles.pdf', help='PDF filename')
parser.add_argument('--xrange', type=str2pair, default=None, help='range of values along x-axis')
parser.add_argument('--yrange', type=str2pair, default=None, help='range of values along y-axis')
parser.add_argument('--pdf', action='store_true', help='save a PDF file instead')
args = parser.parse_args()

class Style:
    thin = { 'linestyle': '-', 'linewidth': 0.5, 'color': 'k' }
    dashed = { 'linestyle': '--', 'linewidth': 0.5, 'color': 'k' }
    point = { 'color': 'black', 'marker': '.', 'linestyle': 'None', 'zorder': 10 }
    annotate = { 'xytext': (0,5), 'textcoords': 'offset points', 'ha': 'center' }

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

_xi = lambda x: 1 - exp(-x)
_rho = lambda r,x: r*exp(-x)
_norm = lambda r,x: sqrt(r**2 + x**2)
_contour_eq = lambda model: lambda r,x: func[model](r,x) - args.Tratio
_depth_eq = lambda model: lambda x: (func[model](*x) - args.Tratio, ddx[model](*x))
_mesh = lambda a, b, va, vb: (erf(np.linspace(-va, vb, args.N)) + 1)*(b-a)/2 + a
_latex = lambda text: r'$\mathrm{' + text.replace(' ', r'\ ') + r'}$'

func = {
    'Levin': lambda r,x: (1-_xi(x))*(1/_norm(_xi(x), _rho(r,x)) - 1/_norm(2-_xi(x), _rho(r,x))),
    'Rosenthal': lambda r,x: exp(-(x+_norm(r,x))/2)/_norm(r,x)
}

ddx = {
    'Levin': lambda r,x: (1-_xi(x))*((2-_xi(x))/_norm(2-_xi(x), _rho(r,x))**3 - _xi(x)/_norm(_xi(x), _rho(r,x))**3),
    'Rosenthal': lambda r,x: -exp(-(x+_norm(r,x))/2)*(x/_norm(r,x)**3 + (1 + x/_norm(r,x))/2/_norm(r,x))
}

if args.debug:
    X = np.linspace(-3, 3, args.N)
    for m in func.keys():
        plt.plot(X, func[m](1,X), label='func')
        plt.plot(X, ddx[m](1,X), label='der')
        plt.title(m)
        plt.legend()
        plt.show()

logN = np.log10(args.N)
small, large = 10*np.finfo(float).eps, 1e2*log(2/args.Tratio)**(1/3)
depthAtZeroX = (args.Tratio)**(-1/3)
lengthApprox = 1.5*log(2/args.Tratio)
if args.verbose:
    print(f'First estimation: length = {lengthApprox:.5g}, depth = {depthAtZeroX:.5g}')

for label in [ 'Levin 2008', 'Rosenthal 1946' ]:
    model = label.split(' ')[0]
    xmin = root_scalar(partial(_contour_eq(model), 0), bracket=[-large, -small]).root
    xmax = root_scalar(partial(_contour_eq(model), 0), bracket=[small, large]).root
    length = xmax - xmin
    depth, xdepth = root(_depth_eq(model), [depthAtZeroX, 0], method='lm').x
    if args.verbose:
        print(f'{label}: {xmin:.5g} < x < {xmax:.5g}, length = {length:.5g}, depth = {depth:.5g}')
    X = _mesh(xmin, xmax, logN+1, logN+1.5)
    Y = [ root_scalar(_contour_eq(model), args=x, bracket=[small, large]).root for x in X ]
    plt.plot(X, -np.array(Y), label=_latex(label))
    plt.plot([xmin, xmax], [0, 0], **Style.point)
    plt.plot(xdepth, -depth, **Style.point)
    if label.startswith('R'):
        plt.annotate('$\hat{x}_\mathrm{min}$', (xmin, 0), **Style.annotate)
        plt.annotate('$\hat{x}_\mathrm{max}$', (xmax, 0), **Style.annotate)
        plt.plot([xdepth, 0], [-depth, -depth], **Style.dashed)
        plt.plot(0, -depth, **Style.point)
        plt.annotate('$-\hat{d}$', (0, -depth), **Style.annotate)

plt.axhline(y=0, **Style.thin)
plt.axvline(x=0, **Style.thin)
plt.axis('scaled')
plt.ylim(top=depth/4)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{z}$', rotation=0)
plt.legend()
if args.xrange:
    plt.xlim(args.xrange)
if args.yrange:
    plt.ylim(args.yrange)

plt.tight_layout(pad=1)
if args.pdf:
    if args.verbose:
        print(f' -- Save to {args.output}')
    plt.savefig(args.output, bbox_inches='tight')
else:
    plt.show()

