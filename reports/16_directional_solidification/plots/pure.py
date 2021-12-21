#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root
from scipy.special import erf
from numpy import pi, sqrt, infty, fabs
from functools import partial

parser = argparse.ArgumentParser(
    description='Solver for analyzing the Mullins--Sekerka instability of a pure substance')

modes = { 'f': 'fixedV', '2': 'kV', 'n': None }
class ParseMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, modes[values[0]])

str2pair = lambda s: [float(item) for item in s.split(':')]

parser.add_argument('mode', choices=[*modes.keys(), *modes.values()], action=ParseMode, help='Execution mode')
parser.add_argument('-N', type=int, default=100, help='number of points')
parser.add_argument('-V', type=float, default=0.02, help='capillary length/diffusion length')
parser.add_argument('-s', '--figsize', type=str2pair, default='5:4', help='figure size')
parser.add_argument('-a', '--asymptotics', action='store_true', help='plot the asymptotics as well')
parser.add_argument('-l', '--log', action='store_true', help='use log scale for V')
parser.add_argument('-w', '--wavelength', action='store_true', help='use wavelength instead of wavenumber')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
parser.add_argument('--krange', type=str2pair, default=None, help='range of wavenumbers')
parser.add_argument('-o', '--output', type=str, default=None, help='PDF filename')
parser.add_argument('--pdf', action='store_true', help='save a PDF file instead')
args = parser.parse_args()

class Style:
    dotted = { 'linestyle': ':', 'linewidth': 0.5, 'color': 'k' }
    dashed = { 'linestyle': '--', 'linewidth': 0.5, 'color': 'k' }
    thin = { 'linestyle': '-', 'linewidth': 0.5, 'color': 'k' }
    thick = { 'linestyle': '-', 'linewidth': 1, 'color': 'k', 'zorder': 10 }
    unstable = { 'label': r'$\mathrm{unstable\ region}$', 'hatch': 'XX', 'color':'none', 'edgecolor': 'gray' }
    gray = { 'color': 'lightgray' }
    point = { 'color': 'black', 'marker': '.', 'linestyle': 'None', 'zorder': 10 }
    annotate = { 'xytext': (0,3), 'textcoords': 'offset points' }

### Global constants
almost_one = 1 - 10*np.finfo(float).eps
factorY = 1.1
Vmax = 1

### Parameter-dependent constants
logN = np.log10(args.N)

### Auxiliary functions in the a_0(k) relation
_L = lambda a,k: 1/2 + sqrt(1/4 + a + k**2)
_L_a = lambda a,k: 1/(2*_L(a,k) - 1)
_L_kk = lambda a,k: 1/(2*_L(a,k) - 1)
_l = lambda k: _L(0,k)
_S = lambda a,k: -1/2 + sqrt(1/4 + a + k**2)
_S_a = lambda a,k: 1/(2*_S(a,k) + 1)
_S_kk = lambda a,k: 1/(2*_S(a,k) + 1)
_s = lambda k: _S(0,k)

### Estimations for k
_kmax = lambda v: sqrt((3**-0.5 - 1/2)*(1-v))   # for V->1
_k0 = lambda v: (2*v)**-0.5                     # for V->0

### Formulas for marginal stability (a_0 = 0)
_k_0_eq = lambda k,v: 1 - _l(k) + v*k**2*(_l(k) + _s(k))

### Formulas for finding the most unstable wavenumber
_k_max_eq = lambda k,v: (1+v-6*v*k**2)/2/v/(1-2*v*k**2) - sqrt((1-2*v*k**2)/v)
_a_max = lambda k,v: (1-2*v*k**2)/4/v - k**2 - 1/4

### Other functions
interval_mesh = lambda a, b, va, vb: (erf(np.linspace(-va, vb, args.N)) + 1)*(b-a)/2 + a
make_log = lambda f: np.log10(f) if args.log else f
k2k = lambda k,v: 2*pi/k/v if args.wavelength else k
klabel = lambda p='', s='': f'${p}'+r'\hat\lambda'+f'{s}$' if args.wavelength else f'${p}k{s}$'
add_tmargin = lambda a,b: (a, b+(factorY-1)*(b-a))

np.seterr(all='raise')  # Consider warnings as errors

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

if args.verbose:
    print(f'Vmax = {Vmax:.5g} at k -> 0')

### Mode 1: Amplification rate (a_0) vs wave number (k)
if args.mode == modes['f']:
    if args.verbose:
        print(f' -- Plot a0(k) for given V = {args.V}')

    omax = 1/8/args.V
    _k2k = lambda k: k2k(k, args.V)

    # Estimation for the neutral stability wavenumber
    k0 = _k0(args.V)
    if args.verbose:
        print(f'k_0 = {k0:.5g}')

    # Equation for a(k) and its derivatives w.r.t. a_0 and k^2
    _a_eq = lambda a,k: a + 1 - _L(a,k) + args.V*k**2*(_L(a,k) + _S(a,k))
    _a_eq_a = lambda a,k: 1 - _L_a(a,k) + args.V*k**2*(_L_a(a,k) + _S_a(a,k))
    _amin = lambda k: almost_one*(-k**2 - 1/4)
    _amean = lambda k: -k**2*(1+args.V) + (args.V*k**2)**2

    if args.krange:
        K = np.geomspace(*args.krange, args.N)
    else:
        K = np.geomspace(_kmax(args.V)/4, 1.5*k0, args.N)

    m0 = _a_eq_a(_amin(K), K) < 0   # if amin < amean, where a2 < amean < a1
    m2 = K < 1/2                    # if a2 is real
    A1, Amean, A2 = [], [], []
    for i, k in enumerate(K):
        if m0[i]:
            Amean.append(_amean(k))
            A1.append(root_scalar(_a_eq, args=k, bracket=[Amean[-1], omax]).root)
            if m2[i]:
                A2.append(root_scalar(_a_eq, args=k, bracket=[_amin(k), Amean[-1]]).root)
        else:
            A1.append(root_scalar(_a_eq, args=k, bracket=[_amin(k), omax]).root)
    A1, Amean, A2 = np.array(A1), np.array(Amean), np.array(A2)
    amax = np.max(A1)

    ### 1. Plot a1(k)
    plt.plot(_k2k(K), A1)
    plt.semilogx()
    plt.xlabel(klabel())
    plt.ylabel(r'$a_0$', rotation=0)
    plt.axhline(y=0, **Style.thin)

    ### 2. Plot a2(k)
    if args.verbose:
        plt.plot(_k2k(K[m2]), A2)
        plt.plot(_k2k(K[m0]), Amean, **Style.dotted)
        plt.fill_between(_k2k(K[m0]), _amin(K)[m0], Amean, **Style.gray)
        plt.axvline(_k2k(1/2), **Style.thin)

    ### 3. Draw k_m and k_0
    k_max = root_scalar(_k_max_eq, args=args.V, bracket=[_kmax(args.V)/2, almost_one*k0]).root
    a_max = _a_max(k_max, args.V)
    k_0 = root_scalar(_k_0_eq, args=args.V, bracket=[k_max, k0]).root

    if args.verbose:
        plt.plot(_k2k(k_max), a_max, **Style.point)
        plt.annotate(klabel(s='_m'), (_k2k(k_max), a_max), **Style.annotate)
        print(f'Maximum point: k = {k_max:.5g}, a = {a_max:.5g}')
        plt.plot(_k2k(k_0), 0, **Style.point)
        plt.annotate(klabel(s='_0'), (_k2k(k_0), 0), **Style.annotate)
        print(f'Neutral stability point: k = {k_0:.5g}')
        plt.ylim(add_tmargin(np.min(np.r_[A1,A2]), np.max(A1)))

    if args.asymptotics:
        _a1_small = lambda k: -sqrt(1-args.V)*k - (args.V+1)*k**2/2
        _a2_small = lambda k: sqrt(1-args.V)*k - (args.V+1)*k**2/2
        plt.plot(_k2k(K[m2]), _a1_small(K[m2]), **Style.dashed)
        plt.plot(_k2k(K[m2]), _a2_small(K[m2]), **Style.dashed)
        if args.V > 0.5:
            _a1_large = lambda k: -k**2 - 1/4 + (2*args.V)**-2
        else:
            _a1_large = lambda k: -1/2 + (1-2*args.V*k**2)/2*sqrt(2/args.V-2)
        plt.plot(_k2k(K[~m2]), _a1_large(K[~m2]), **Style.dashed)

### Mode 2: Stability diagram in the (V,k) coordinates
elif args.mode == modes['2']:
    if args.verbose:
        print(f' -- Plot the (V,k) stability diagram')

    if (args.V <= 0 or args.V >= Vmax):
        print('The planar front is unconditionally stable.')
        sys.exit()

    V = interval_mesh(0, Vmax, 1.8, logN+1)

    K_most = np.array([
        root_scalar(_k_max_eq, args=v, bracket=[_kmax(v)/2, almost_one*_k0(v)]).root for v in V ])
    A_most = _a_max(K_most, V)
    K_0 = np.array([ root_scalar(_k_0_eq, args=v, bracket=[k, _k0(v)]).root for v,k in zip(V,K_most) ])
    K_min = 2e-2*np.ones_like(K_0)

    if args.wavelength:
        K_min = K_min/V
        plt.ylim(np.min(k2k(K_0,V))/factorY, k2k(K_min,V)[0])
    else:
        plt.ylim(K_min[0], factorY*K_0[0])

    plt.plot(V, k2k(K_0,V))
    plt.plot(V, k2k(K_most,V), label=r'$\mathrm{the\ most\ unstable}$')
    plt.fill_between(V, k2k(K_min,V), k2k(K_0,V), **Style.unstable)
    plt.loglog() if args.log else plt.semilogy()
    plt.xlabel(r'$\hat{V}$')
    plt.ylabel(klabel(), rotation=0)
    plt.xlim(V[0], Vmax)
    plt.legend()

    if args.asymptotics:
        V1, V2 = np.split(V, 2)
        plt.plot(V1, k2k((6*V1)**-0.5, V1), **Style.dashed)
        plt.plot(V1, k2k(_k0(V1), V1), **Style.dashed)
        plt.plot(V2, k2k(sqrt((1-V2)/3), V2), **Style.dashed)
        plt.plot(V2, k2k(_kmax(V2), V2), **Style.dashed)

if args.mode:
    filename = args.output if args.output else f'{args.mode}.pdf'
    plt.tight_layout()
    if args.pdf:
        if args.verbose:
            print(f' -- Save to {filename}')
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()

