#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root
from scipy.special import erf
from numpy import sqrt, exp, log, arctan2, cos, tan
from functools import partial

parser = argparse.ArgumentParser(description='Draw a solid--liquid interface of a melt pool')

modes = { 'b': 'boundary', 'g': 'gradient', 's': 'speed', 'r': 'rate' }
class ParseMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, modes[values[0]])

str2pair = lambda s: [float(item) for item in s.split(':')]

parser.add_argument('mode', choices=[*modes.keys(), *modes.values()], action=ParseMode, help='Execution mode')
parser.add_argument('-N', type=int, default=100, help='number of points along each axis')
parser.add_argument('-T', '--Tratio', type=float, default=0.05, help='T_M/T_0')
parser.add_argument('-s', '--figsize', type=str2pair, default='5:4', help='figure size')
parser.add_argument('-a', '--arc', action='store_true', help='use the arc-length coordinate instead of x')
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
_gradT = lambda model: lambda r,x: _norm(ddx[model](r,x), ddr[model](r,x))
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

ddr = {
    'Levin': lambda r,x: _rho(r,x)*(1-_xi(x))**2*(1/_norm(2-_xi(x), _rho(r,x))**3 - 1/_norm(_xi(x), _rho(r,x))**3),
    'Rosenthal': lambda r,x: -exp(-(x+_norm(r,x))/2)*r/_norm(r,x)**2*(1/_norm(r,x) + 1/2)
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
if args.verbose and args.mode == modes['b']:
    print(f'First estimation: length = {lengthApprox:.5g}, depth = {depthAtZeroX:.5g}')

for n, label in enumerate([ 'Levin 2008', 'Rosenthal 1946' ]):
    model = label.split(' ')[0]

    # 1. Find the melt-pool boundary
    xmin = root_scalar(partial(_contour_eq(model), 0), bracket=[-large, -small]).root
    xmax = root_scalar(partial(_contour_eq(model), 0), bracket=[small, large]).root
    length = xmax - xmin
    xminmax = np.array([xmin, xmax])
    X = _mesh(xmin, xmax, logN+1, logN+1.5)
    R = np.array([ root_scalar(_contour_eq(model), args=x, bracket=[small, large]).root for x in X ])

    # 2. Find the deepest point
    depth, xdepth = root(_depth_eq(model), [depthAtZeroX, 0], method='lm').x
    sol, fus = X <= xdepth, X > xdepth

    # 3. Find temperature gradient
    gradT_depth = _gradT(model)(depth, xdepth)
    gradT_minmax = _gradT(model)(np.zeros(2), xminmax)
    GradT = _gradT(model)(R,X)

    # 4. Find the normal speed of the melt-pool boundary
    Phi = arctan2(ddr[model](R,X), ddx[model](R,X))
    CosPhi = cos(arctan2(ddr[model](R,X), ddx[model](R,X)))
    rate_minmax = np.array([1,-1])*gradT_minmax

    # 5. Find the coordinate along the curve using the midpoint integration rule
    DX = np.zeros_like(X)
    DX[0], DX[-1] = X[0] - xmin, xmax - X[-1]
    DX[1:-2] = (X[2:-1] - X[0:-3])/2
    S = np.cumsum(sqrt(1 + 1/tan(Phi)**2)*DX)
    arc_length = S[-1]
    S /= arc_length

    ### Mode 1: Melt-pool boundary
    if args.mode == modes['b']:
        if args.verbose:
            print(f'{label:>14s}: {xmin:.5g} < x < {xmax:.5g}, length = {length:.5g}, '
                + f'depth = {depth:.5g}, boundary length = {arc_length:.5}')
        plt.plot(X, -R, label=_latex(label))
        plt.plot(xminmax, np.zeros(2), **Style.point)
        plt.plot(xdepth, -depth, **Style.point)
        plt.axvline(**Style.thin)
        plt.axis('scaled')
        plt.ylim(top=depth/4)
        plt.ylabel(r'$\hat{z}$', rotation=0)
        if label.startswith('R'):
            plt.annotate('$\hat{x}_\mathrm{min}$', (xmin, 0), **Style.annotate)
            plt.annotate('$\hat{x}_\mathrm{max}$', (xmax, 0), **Style.annotate)
            plt.plot([xdepth, 0], [-depth, -depth], **Style.dashed)
            plt.plot(0, -depth, **Style.point)
            plt.annotate('$-\hat{d}$', (0, -depth), **Style.annotate)

    ### Mode 2: Temperature gradient along the melt-pool boundary
    elif args.mode == modes['g']:
        GradT1, GradT2 = GradT[sol], GradT[fus]
        if args.verbose:
            print(f'{label:>14s}: temperature gradient at solidification: ' +
                f'min = {np.min(GradT1):.5g}, max = {np.max(GradT1):.5g} at fusion: ' +
                f'min = {np.min(GradT2):.5g}, max = {np.max(GradT2):.5g}')
        if args.arc:
            X1, X2 = S[sol], S[fus]
            xminmax = [0,1]
            xdepth = np.interp(xdepth, X, S)
            for x in [0,1]:
                plt.axvline(x=x, **Style.thin)
        else:
            X1, X2 = X[sol], X[fus]
        plt.plot(xdepth, gradT_depth, **Style.point)
        X1, GradT1 = np.append(X1, xdepth), np.append(GradT1, gradT_depth)
        plt.plot(X1, GradT1, label=_latex(label))
        plt.plot(X2, GradT2, color=f'C{n}', linestyle='--')
        plt.plot(xminmax, gradT_minmax, **Style.point)
        plt.ylabel(r'$\hat{\nabla}\hat{T}$', rotation=0, ha='right')

    ### Mode 3: Speed of the melt-pool boundary (positive is for solidification)
    elif args.mode == modes['s']:
        for y in [-1,1]:
            plt.axhline(y=y, **Style.thin)
        if args.arc:
            xminmax = [0,1]
            xdepth = np.interp(xdepth, X, S)
            X = S
        plt.plot(X, CosPhi, label=_latex(label))
        plt.plot(xdepth, 0, **Style.point)
        plt.plot(xminmax, [1,-1], **Style.point)
        plt.ylabel(r'$\cos\phi$', rotation=0, ha='right')

    ### Mode 4: Cooling/heating rate along the melt-pool boundary
    elif args.mode == modes['r']:
        if args.verbose:
            print(f'{label:>14s}: cooling rate: max = {rate_minmax[0]:.5g}' +
                f' heating rate: max = {-rate_minmax[1]:.5g}')
        if args.arc:
            xminmax = [0,1]
            xdepth = np.interp(xdepth, X, S)
            X = S
        plt.plot(X, CosPhi*GradT, label=_latex(label))
        plt.plot(xdepth, 0, **Style.point)
        plt.plot(xminmax, rate_minmax, **Style.point)
        plt.ylabel(r'$\hat{\nabla}\hat{T}\cos\phi$', rotation=0, ha='right')

plt.axhline(**Style.thin)
plt.legend()
if args.arc and args.mode != modes['b']:
    plt.xlabel(r'$s$')
else:
    plt.xlabel(r'$\hat{x}$')

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

