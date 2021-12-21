#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root
from scipy.special import erf
from numpy import pi, sqrt, infty, fabs
from functools import partial
from termcolor import colored

parser = argparse.ArgumentParser(
    description='Solver for analyzing the Mullins--Sekerka instability of a binary mixture')

modes = { 'b': 'bifurcation', 'f': 'fixedGV', '2': '2d', '3': '3d', 'n': None }
class ParseMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, modes[values[0]])

str2pair = lambda s: [float(item) for item in s.split(':')]

parser.add_argument('mode', choices=[*modes.keys(), *modes.values()], action=ParseMode, help='Execution mode')
parser.add_argument('-N', type=int, default=100, help='number of points')
parser.add_argument('-k', '--kratio', type=float, default=1, help='k_S/k_L')
parser.add_argument('-D', '--Dratio', type=float, default=0, help='D_S/D_L')
parser.add_argument('-K', type=float, default=0.5, help='partition coefficient')
parser.add_argument('-V', type=float, default=0.02, help='capillary length/diffusion length')
parser.add_argument('-G', type=float, default=0.01, help='capillary length/thermal length')
parser.add_argument('-s', '--size', type=float, default=5, help='figure size')
parser.add_argument('-a', '--asymptotics', action='store_true', help='plot the asymptotics as well')
parser.add_argument('-l', '--log', action='store_true', help='use log scale for V and G')
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
    surface = { 'rstride': 1, 'cstride': 1, 'color': 'gray', 'edgecolor': 'none', 'alpha': 0.5 }
    annotate = { 'xytext': (0,3), 'textcoords': 'offset points' }

### Global constants
almost_one = 1 - 10*np.finfo(float).eps
factorY = 1.1
kmin, kmax = 1e-5, 1e5  #TODO: provide some reasonable estimates

### Parameter-dependent constants
Vmax = fabs(1-args.K)/args.K
etaK = args.Dratio*args.K
logN = np.log10(args.N)

### Auxiliary functions in the a_0(k) relation
_L = lambda a,k: 1/2 + sqrt(1/4 + a + k**2)
_L_a = lambda a,k: 1/(2*_L(a,k) - 1)
_L_kk = lambda a,k: 1/(2*_L(a,k) - 1)
_l = lambda k: _L(0,k)
_S = lambda a,k: -1/2 + sqrt(1/4 + args.Dratio*a + (args.Dratio*k)**2)
_S_a = lambda a,k: args.Dratio/(2*_S(a,k) + 1)
_S_kk = lambda a,k: args.Dratio**2/(2*_S(a,k) + 1)
_s = lambda k: _S(0,k)

### Formulas for calculating beta, V, calG, G
_beta = lambda v: v/fabs(1-args.K)
_V = lambda beta: beta*fabs(1-args.K)
_calG = lambda g,v: 1 - 2*g/(args.kratio+1)/v/fabs(1-args.K)
_G = lambda calG,v: (args.kratio+1)*v*fabs(1-args.K)*(1 - calG)/2

### Formulas for marginal stability (a_0 = 0)
_denom = lambda k: (_l(k)-1) + args.K*(1+_s(k))
_calG_kv = lambda k,v: _beta(v)*k**2 + args.K*(1+_s(k))/_denom(k)
_G_kv = lambda k,v: _G(_calG_kv(k,v), v)

### Formulas for finding bifurcation points
_beta_bif = lambda k: args.K/_denom(k)**2*(
    (1+_s(k))/(2*_l(k)-1) - args.Dratio**2*(_l(k)-1)/(2*_s(k)+1) )
_V_bif = lambda k: _V(_beta_bif(k))

### The 2D equation for finding the most unstable wavenumber
__f = lambda a,k,g,v: a + args.K*(1+_S(a,k)) - (_calG(g,v) - _beta(v)*k**2)*(
    _L(a,k)-1 + args.K*(1+_S(a,k)) )
__f_kk = lambda a,k,g,v: args.K*_S_kk(a,k) - (_calG(g,v) - _beta(v)*k**2)*(
    _L_kk(a,k) + args.K*_S_kk(a,k)) + _beta(v)*(_L(a,k)-1 + args.K*(1+_S(a,k)))
_most_eq = lambda a,k,g,v: (__f(a,k,g,v), __f_kk(a,k,g,v))

### Other functions
interval_mesh = lambda a, b, va, vb: (erf(np.linspace(-va, vb, args.N)) + 1)*(b-a)/2 + a
make_log = lambda f: np.log10(f) if args.log else f
k2k = lambda k,v: 2*pi/k/v if args.wavelength else k
_k2k = lambda k: k2k(k, args.V)
klabel = lambda p='', s='': f'${p}'+r'\hat\lambda'+f'{s}$' if args.wavelength else f'${p}k{s}$'
add_tmargin = lambda a,b: (a, b+(factorY-1)*(b-a))

def error(msg):
    print(colored(msg, 'red'), file=sys.stderr)

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
        'figure.figsize': [5,4]
    }
    plt.rcParams.update(params)

### Calculate Gmax
_k_star_eq = lambda k: _calG_kv(k, _V_bif(k)) + _beta_bif(k)*k**2 - 1
k_star = root_scalar(_k_star_eq, bracket=[1e-3, 1e3]).root
v_star = _V(_beta_bif(k_star))
Gmax = _G_kv(k_star, v_star)
if args.verbose:
    print(f'Vmax = {Vmax:.5g} at G = 0 and k -> 0')
    print(f'Gmax = {Gmax:.5g} at V = {v_star:.5g} and k = {k_star:.5g}')

### Mode 1: Bifurcation diagram in the (V,G) coordinates
### Stable region corresponds to the MS stability for any wave numbers
if args.mode == modes['b']:
    factor = 1 if args.Dratio > 0 else 2    # NB: asymptotics for Dratio=0 differ!

    K = np.logspace(-1, 1, args.N)*k_star
    V = _V_bif(K)
    G = _G_kv(K, V)

    if args.verbose:
        fig, axs = plt.subplots(ncols=2, figsize=args.size*np.array((2.2,1)))
        axs[0].plot(k2k(K,V), V)
        axs[0].plot(k2k(k_star,v_star), v_star, **Style.point)
        axs[0].set_xlabel(klabel(s='_\mathrm{bif}'))
        axs[0].set_ylabel(r'$\hat{V}$', rotation=0, color='C0')

        ax = axs[0].twinx()
        ax.plot(k2k(K,V), G, color='C1')
        ax.plot(k2k(k_star,v_star), Gmax, **Style.point)
        ax.set_ylabel(r'$\hat{G}$', rotation=0, color='C1')

        if args.log:
            axs[0].loglog(); ax.loglog()
        else:
            axs[0].semilogx(); ax.semilogx()


        if args.asymptotics:
            K1, K2 = np.split(K, 2)
            V1, V2 = np.split(V, 2)
            # Asymptotic approximation of v(k) for small v
            _vsmall = lambda k: _V(factor*args.K*(1+args.Dratio)/4/k**3/(1+args.K*args.Dratio)**2)
            axs[0].plot(k2k(K2,V2), _vsmall(K2), **Style.dashed)

    K = np.logspace(-logN+1, logN-1, args.N)*k_star
    V = _V_bif(K)
    G = _G_kv(K,V)

    ax = axs[1] if args.verbose else plt
    ax.plot(V, G)
    ax.fill_between(V, 0*G, G, **Style.unstable)
    ax.legend()
    ax.plot(v_star, Gmax, **Style.point)
    ax.annotate(r'$\hat{G}_\mathrm{max}$', (v_star, Gmax), **Style.annotate)

    if args.verbose:
        axs[1].set_xlabel(r'$\hat{V}$')
        axs[1].set_ylabel(r'$\hat{G}$', rotation=0)
    else:
        plt.xlabel(r'$\hat{V}$')
        plt.ylabel(r'$\hat{G}$', rotation=0)

    if args.log:
        plt.loglog()
    else:
        if args.verbose:
            axs[1].set_xlim(0, Vmax)
            axs[1].set_ylim(add_tmargin(0, Gmax))
        else:
            plt.xlim(0, Vmax)
            plt.ylim(add_tmargin(0, Gmax))

    if args.asymptotics:
        _calG_asym1 = lambda k: args.K/(1 + args.K*args.Dratio)*(
            args.Dratio + (1 + (1+args.Dratio)/2/(1+args.K*args.Dratio))*factor/2/k )
        _calG_asym2 = lambda k: 1 - 1*k**4 # TODO: calculate the coefficient
        K1, K2 = np.split(K, 2)

        V1 = _vsmall(K2)
        G1 = _G(_calG_asym1(K2), V1)
        if args.verbose:
            axs[1].plot(V1, G1, **Style.dashed)
        else:
            plt.plot(V1, G1, **Style.dashed)

### Mode 2: Amplification rate (a_0) vs wave number (k)
elif args.mode == modes['f']:
    beta, calG = _beta(args.V), _calG(args.G, args.V)
    omax = 10*(1 + sqrt(args.K*etaK))**2

    # Simple estimate for the reference wavenumber (depends on G, V, kratio, Dratio, K)
    _k0 = lambda g,v: sqrt(max(0, _calG(g,v)*(1+etaK) - etaK)/_beta(v)/(1+etaK))
    k0 = _k0(args.G, args.V)
    if args.verbose:
        print(f'k_0 = {k0:.5g}')

    # Equation for a_0(k) and its derivatives w.r.t. a_0 and k^2
    _a_eq = lambda a,k: a + args.K*(1+_S(a,k)) - (calG - beta*k**2)*(
        _L(a,k)-1 + args.K*(1+_S(a,k)) )
    _a_eq_a = lambda a,k: 1 + args.K*_S_a(a,k) - (calG - beta*k**2)*(
        _L_a(a,k) + args.K*_S_a(a,k) )
    _a_eq_kk = lambda a,k: args.K*_S_kk(a,k) - (calG - beta*k**2)*(
        _L_kk(a,k) + args.K*_S_kk(a,k)) + beta*(_L(a,k)-1 + args.K*(1+_S(a,k)))
    # Lower boundary for a_0
    _aminL = lambda k: -k**2 - 1/4
    _aminS = lambda k: -args.Dratio*k**2 - 1/4/args.Dratio if args.Dratio > 0 else -infty
    _amin = lambda k: almost_one*np.maximum(_aminL(k), _aminS(k))
    # For Dratio = 0: _amean = lambda k: -k**2 - (1 - (calG-beta*k**2)**2)/4
    _amean = lambda k: root_scalar(_a_eq_a, args=k, bracket=[_amin(k), omax]).root

    ### 1. Find a critical point (where a1(k)=a2(k)); a_0 is complex behind this point
    if args.G > 0:
        # Try the exact solution for Dratio=0 as initial guess
        A = calG + 2*args.K - 1
        D = (A-2/beta)**2 + (2*args.K)**2 - A**2
        if D > 0:
            B = A-2/beta + sqrt(D)
            if B > 0:
                k_crit = sqrt(B/beta)
                _a_crit = lambda k: -k**2 - (1 - (calG-beta*k**2)**2)/4
                a_crit = _a_crit(k_crit)
        try:
            k_crit
        except NameError:
            error('Failed to use a good initial guess!')
            a_crit, k_crit = 1, 1
        if args.verbose:
            print(f'Initial guess: k = {k_crit:.5g}, a = {a_crit:.5g}')
        try:
            _ak_crit_eq = lambda x: (_a_eq(*x), _a_eq_a(*x))
            res = root(_ak_crit_eq, [a_crit, k_crit])
            if not res.success:
                raise ValueError(res.message)
            a_crit, k_crit = res.x
            if args.verbose:
                plt.plot(_k2k(k_crit), a_crit, **Style.point)
                shift = Style.annotate['xytext'][1]
                style = Style.annotate | { 'xytext': (shift, shift) }
                plt.annotate(klabel(s='_c'), (_k2k(k_crit), a_crit), **style)
                print(f'Critical point: k = {k_crit:.5g}, a = {a_crit:.5g}')
        except (ValueError, FloatingPointError) as err:
            error(f'Failed to find a critical point: {err}')
    else:
        k_crit = 0

    ### 2. Find the extremum point and roots of a1(k)
    try:
        _ak_max_eq = lambda x: (_a_eq(*x), _a_eq_kk(*x))
        a0 = -_a_eq(0, k0)  # Initial guess in the quasi-stationary approximation
        # NB: Levenberg-Marquardt algorithm is used here as more stable
        res = root(_ak_max_eq, [a0, k0], method='lm')
        if not res.success:
            raise ValueError(res.message)
        a_max, k_max = res.x
        k_max = np.abs(k_max)   # Since a negative solution is equivalent to a positive one
        if args.verbose:
            plt.plot(_k2k(k_max), a_max, **Style.point)
            plt.annotate(klabel(s='_m'), (_k2k(k_max), a_max), **Style.annotate)
            print(f'Maximum point: k = {k_max:.5g}, a = {a_max:.5g}')
        Kzero = []
        if a_max > 0:
            Kzero.append(root_scalar(partial(_a_eq, 0), x0=k_max, x1=k_max*2).root)
            if args.G > 0:
                Kzero.append(root_scalar(partial(_a_eq, 0), x0=k_max, x1=k_max/2).root)
            if args.verbose:
                plt.plot(_k2k(Kzero), np.zeros_like(Kzero), **Style.point)
                I = [2,1] if args.G > 0 else [0]
                for i,k in zip(I,Kzero):
                    shift = Style.annotate['xytext'][1]
                    style = Style.annotate | { 'xytext': (-4*shift if i==1 else 0, shift) }
                    plt.annotate(klabel(s=f'_{i}'), (_k2k(k), 0), **style)
                print(f'Neutral stability points: k =', ', '.join(f'{k:.5g}' for k in Kzero))

        else:
            print('The planar front is unconditionally stable.')
    except (ValueError, FloatingPointError) as err:
        error(f'Failed to find a maximum point: {err}')

    ### 3. Create a mesh and find solutions on it
    K = np.geomspace(*args.krange, args.N) if args.krange else np.logspace(-1.5, 0.5, args.N)*k0
    if k_crit > K[0]:
        # Refine the mesh near the critical point
        K = np.r_[K[K<k_crit], interval_mesh(k_crit, K[-1], logN+1, 0)]
    m0 = _a_eq_a(_amin(K), K) < 0           # if amin < amean, where a2 < amean < a1
    m1 = np.zeros_like(K, dtype=bool)       # if a1 (larger) is real
    m2 = np.copy(m1)                        # if a2 (smaller) is real
    A1, Amean, A2 = [], [], []
    for i, k in enumerate(K):
        if m0[i]:
            Amean.append(_amean(k))
            if _a_eq(Amean[-1], k) <= 0:
                m1[i] = True
                A1.append(root_scalar(_a_eq, args=k, bracket=[Amean[-1], omax]).root)
                if _a_eq(_amin(k), k) >= 0:
                    m2[i] = True
                    A2.append(root_scalar(_a_eq, args=k, bracket=[_amin(k), Amean[-1]]).root)
        else:
            if _a_eq(_amin(k), k) <= 0:
                m1[i] = True
                A1.append(root_scalar(_a_eq, args=k, bracket=[_amin(k), omax]).root)
    A1, Amean, A2 = np.array(A1), np.array(Amean), np.array(A2)
    amax = np.max(A1)

    ### 4. Plot a1(k)
    plt.plot(_k2k(K[m1]), A1)
    plt.semilogx()
    plt.xlabel(klabel())
    plt.ylabel(r'$a_0$', rotation=0)
    plt.axhline(y=0, **Style.thin)

    ### 5. Plot a2(k)
    if args.verbose:
        plt.plot(_k2k(K[m2]), A2)
        plt.plot(_k2k(K[m0]), Amean, **Style.dotted)
        plt.fill_between(_k2k(K[m0]), _amin(K)[m0], Amean, **Style.gray)
        plt.ylim(add_tmargin(np.min(np.r_[A1,A2]), np.max(A1)))

### Mode 3a: Stability diagram in the (G,k) and (V,k) coordinates
elif args.mode == modes['2']:
    if (args.G < 0 or args.G > Gmax or args.V <= 0 or args.V >= Vmax):
        print('The planar front is unconditionally stable for given G and V.')
        sys.exit()

    fig, axs = plt.subplots(ncols=2, figsize=args.size*np.array((2,1)))

    ### Subfigure A
    if args.verbose:
        print(f' == (V,k) diagram for given G = {args.G} ==')

    ### A1. Find two bifurcation points and create a mesh between them
    if args.G > 0:
        _k_eq = lambda k: _calG_kv(k, _V_bif(k)) - _calG(args.G, _V_bif(k))

        k1 = root_scalar(_k_eq, bracket=[kmin, k_star]).root
        k2 = root_scalar(_k_eq, bracket=[k_star, kmax]).root
        K_bif = np.array([k2, k1])
        V_bif = _V_bif(K_bif)
        V = interval_mesh(*V_bif, logN+2, logN+1)
    else:
        K_bif = [np.inf, 0]
        V_bif = [0, Vmax]
        V = interval_mesh(*V_bif, -1, logN)

    if args.verbose:
        for i,s in zip(range(2), ['Min', 'Max']):
            print(f'{s} unstable V = {V_bif[i]:.5g} with k = {K_bif[i]:.5g}')

    ### A2. Find the most unstable curve
    if args.G > 0:
        # NB: it is crucial that _k_guess(v) is a straight line on a log-log plot
        _k_guess = lambda v: np.exp(np.interp(np.log(v), np.log(V_bif), np.log(K_bif)))
    else:
        # NB: kc is not the best guess since kc~V^1/3, but k_most~V^1/2
        _k_guess = np.vectorize(lambda v: root_scalar(lambda k: _V_bif(k)-v, bracket=[kmin, kmax]).root)

    A_most, K_most = np.array([
        root(lambda x: _most_eq(*x, args.G, v), [0, _k_guess(v)], method='lm').x for v in V ]).T

    ### A3. Find the boundary of the MS instability
    _k_eq = lambda k,v: _calG_kv(k,v) - _calG(args.G, v)
    K2 = np.array([ root_scalar(_k_eq, args=v, bracket=[k, kmax]).root for k,v in zip(K_most,V) ])
    if args.G > 0:
        K1 = np.array([ root_scalar(_k_eq, args=v, bracket=[kmin, k]).root for k,v in zip(K_most,V) ])
    else:
        K1 = 2e-2*np.ones_like(K2)
        if args.wavelength:
            K1 = K1/V
            axs[0].ylim(np.min(k2k(K2,V))/factorY, k2k(K1,V)[0])
        else:
            axs[0].ylim(K1[0], factorY*K2[0])

    ### A4. Plot the stability diagram
    axs[0].plot(V, k2k(K1,V))
    axs[0].plot(V, k2k(K2,V), color='C0')
    axs[0].plot(V, k2k(K_most,V), label=r'$\mathrm{the\ most\ unstable}$')
    axs[0].fill_between(V, k2k(K1,V), k2k(K2,V), **Style.unstable)
    axs[0].loglog() if args.log else axs[0].semilogy()
    axs[0].set_xlabel(r'$\hat{V}$')
    axs[0].set_ylabel(klabel(), rotation=0)
    axs[0].legend()

    if args.verbose:
        if args.G > 0:
            axs[0].plot(V_bif, k2k(K_bif,V_bif), **Style.point)

    ### Subfigure B
    if args.verbose:
        print(f' == (G,k) diagram for given V = {args.V} ==')

    ### B1. Find the bifurcation point and create a mesh up to it
    _k_eq = lambda k: _V_bif(k) - args.V

    k_bif = root_scalar(_k_eq, bracket=[kmin, kmax]).root
    G_bif = _G_kv(k_bif, args.V)
    G = interval_mesh(0, G_bif, 1, logN+1)

    if args.verbose:
        print(f'Max unstable G = {G_bif:.5g} with k = {k_bif:.5g}')

    ### B2. Find the most unstable curve
    A_most, K_most = np.array([
        root(lambda x: _most_eq(*x, g, args.V), [0, k_bif], method='lm').x for g in G ]).T

    ### B3. Find the boundary of the MS instability
    _k_eq = lambda k,g: _calG_kv(k, args.V) - _calG(g, args.V)

    K1 = np.array([ root_scalar(_k_eq, args=g, bracket=[kmin, k]).root for k,g in zip(K_most,G) ])
    K2 = np.array([ root_scalar(_k_eq, args=g, bracket=[k, kmax]).root for k,g in zip(K_most,G) ])

    ### B4. Plot the stability diagram
    axs[1].plot(G, _k2k(K1))
    axs[1].plot(G, _k2k(K2), color='C0')
    axs[1].fill_between(G, _k2k(K1), _k2k(K2), **Style.unstable)
    axs[1].plot(G, _k2k(K_most), label=r'$\mathrm{the\ most\ unstable}$')
    axs[1].loglog() if args.log else axs[1].semilogy()
    axs[1].set_xlabel(r'$\hat{G}$')
    axs[1].set_ylabel(klabel(), rotation=0)
    axs[1].legend()

    if args.verbose:
        axs[1].plot(G_bif, _k2k(k_bif), **Style.point)

### Mode 3b: Stability diagram in the (V,G,k) coordinates
elif args.mode == modes['3']:
    from mpl_toolkits import mplot3d
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    ### 1. Create a 2D mesh in the (G,t) coordinates,
    ###     where t is a normalized distance between the bifurcation points
    G = interval_mesh(0, Gmax, 1, logN+1)
    T = interval_mesh(0, 1, logN+2, logN+1)
    GG, TT = np.meshgrid(G, T)

    ### 2. Find two bifurcation curves and create meshes between them
    _k_eq = lambda k,g: _calG_kv(k, _V_bif(k)) - _calG(g, _V_bif(k))
    K1 = np.array([ root_scalar(_k_eq, args=g, bracket=[kmin, k_star]).root for g in G ])
    K2 = np.array([ root_scalar(_k_eq, args=g, bracket=[k_star, kmax]).root for g in G ])

    K12 = np.array([K2,K1])
    V12 = _V_bif(K12)

    _t2v = lambda t,v1,v2: np.interp(t, [0,1], [v1,v2])
    VV = np.array([ _t2v(T, *v12) for v12 in V12.T ]).T

    ### 3. Find the most unstable surface
    _k_guess = lambda v, v12, k12: np.exp(np.interp(np.log(v), np.log(v12), np.log(k12)))
    KK_guess = np.array([ _k_guess(_t2v(T, *v12), v12, k12) for k12,v12 in zip(K12.T,V12.T) ]).T

    _most_eq_gvk = np.vectorize(lambda g,v,k_guess:
        tuple(root(lambda x: _most_eq(*x,g,v), [0,k_guess], method='lm').x))
    AA_most, KK_most = _most_eq_gvk(GG, VV, KK_guess)

    ### 4. Find the boundary of the MS instability
    _k_eq = lambda k,g,v: _calG_kv(k,v) - _calG(g,v)
    _k1_eq_gvk = np.vectorize(lambda g,v,k_most:
        root_scalar(_k_eq, args=(g,v), bracket=[k_most, kmax]).root)
    _k2_eq_gvk = np.vectorize(lambda g,v,k_most:
        root_scalar(_k_eq, args=(g,v), bracket=[kmin, k_most]).root)

    KK1 = _k1_eq_gvk(GG, VV, KK_most)
    KK2 = _k2_eq_gvk(GG, VV, KK_most)

    ### 5. Set logarithmic scales
    # NB: due to some bug in Matplotlib, ax.set_xscale('log') doesn't work for all axes;
    #   therefore, we have to transform the data manually.
    KK_most = np.log10(k2k(KK_most,VV)); KK1 = np.log10(k2k(KK1,VV)); KK2 = np.log10(k2k(KK2,VV))
    VV = make_log(VV); GG = make_log(GG)
    ax.set_zlabel(klabel(p=r'\log_{10}'))
    if args.log:
        ax.set_xlabel(r'$\log_{10}\hat{V}$'); ax.set_ylabel(r'$\log_{10}\hat{G}$')
    else:
        ax.set_xlabel(r'$\hat{V}$'); ax.set_ylabel(r'$\hat{G}$')

    ### 6. Plot the 3D stability diagram
    ax.plot_surface(VV, GG, KK1, **Style.surface)
    ax.plot_surface(VV, GG, KK2, **Style.surface)
    p = ax.plot_surface(VV, GG, KK_most, cmap='viridis')
    fig.colorbar(p)

    if args.verbose:
        v0 = make_log(v_star); g0 = make_log(Gmax); k0 = np.log10(k2k(k_star,v_star))
        ax.plot3D(v0, g0, k0, **Style.point)
        V1 = make_log(_V_bif(K1)); V2 = make_log(_V_bif(K2)); G = make_log(G)
        K1 = np.log10(k2k(K1,_V_bif(K1))); K2 = np.log10(k2k(K2,_V_bif(K2)))
        ax.plot3D(V1, G, K1, **Style.thick)
        ax.plot3D(V2, G, K2, **Style.thick)

if args.mode:
    filename = args.output if args.output else f'{args.mode}.pdf'
    plt.tight_layout()
    if args.pdf:
        if args.verbose:
            print(f'Saving to {filename}...')
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()
