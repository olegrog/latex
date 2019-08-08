#!/usr/bin/env python3
import sys, argparse, threading, logging, traceback, vgrid
import numpy as np
from functools import partial
from collections import namedtuple

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-U', type=float, default=2e-2, metavar='<float>', help='difference in velocity of plates')
parser.add_argument('-k', '--kn', type=float, default=0.1, metavar='<float>', help='modified Knudsen number')
parser.add_argument('-N1', type=int, default=4, metavar='<int>', help='number of cells in the first domain')
parser.add_argument('-N2', type=int, default=4, metavar='<int>', help='number of cells in the second domain')
parser.add_argument('-g', '--grid', default='polynomial', metavar='uniform|hermite|polynomial|geometric', help='type of the grid along xi_y')
parser.add_argument('-M', type=int, default=8, metavar='<int>', help='number of points along xi_x, xi_z')
parser.add_argument('-My', type=int, default=2, metavar='<int>', help='factor for number of points along xi_y')
parser.add_argument('-q', '--ratio', type=float, default=1.15, metavar='<float>', help='ratio for the geometric grid')
parser.add_argument('--w-min', type=float, default=0.1, metavar='<float>', help='minimum weight for the polynomial grid')
parser.add_argument('--poly', type=float, default=2, metavar='<float>', help='power for the polynomial grid')
parser.add_argument('-e', '--end', type=int, default=int(1e2), metavar='<int>', help='maximum total time')
parser.add_argument('-w', '--width', type=str, default='2.4/np.sqrt(np.pi)', metavar='<expr>', help='width of the right domain in terms of the Knudsen number')
parser.add_argument('-s', '--step', type=float, default=1.0, metavar='<float>', help='reduce timestep in given times')
parser.add_argument('-R', '--radius', type=float, default=4., metavar='<float>', help='radius of the velocity grid')
parser.add_argument('-r', '--refinement', type=float, default=2., metavar='<float>', help='ratio of maximum and minumum cell width')
parser.add_argument('-m1', '--model1', default='lbm', metavar='dvm|lbm', help='type of velocity model in the first domain')
parser.add_argument('-m2', '--model2', default='dvm', metavar='dvm|lbm', help='type of velocity model in the second domain')
parser.add_argument('-l', '--lattice', default='d3q19', metavar='d3q15|19|27', help='type of velocity lattice')
parser.add_argument('-i', '--limiter', default='wide-third', metavar='none|minmod|mc|superbee|...', help='type of limiter')
parser.add_argument('-c', '--correction', default='poly', metavar='none|poly|entropy|lagrange', help='type of correction for construction of discrete Maxwellian')
parser.add_argument('-o', '--reconstruction', default='hermite40', metavar='grad13|hermite40', help='type of reconstruction in the buffer layer')
parser.add_argument('-n', '--norm', default='L2', metavar='L1|L2|Linf', help='norm for distance between solutions')
parser.add_argument('-p', '--plot', type=int, default=10, metavar='<int>', help='plot every <int> steps')
parser.add_argument('--plot-norms', action='store_true', help='plot profiles of the norms instead')
parser.add_argument('--plot-ns', action='store_true', help='plot profiles for Navier--Stokes--Fourier equations')
parser.add_argument('--qflow-factor', type=float, default=40, metavar='<float>', help='plot <float>*qflow instead of qflow')
parser.add_argument('--epsilon', type=float, default=1e-5, metavar='<float>', help='maximum error in test mode')
parser.add_argument('-t', '--tests', action='store_true', help='run some tests instead')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

if args.plot:
    import matplotlib.pyplot as plt

log_level = logging.DEBUG if args.verbose else logging.INFO
logging.basicConfig(level=log_level, format='(%(threadName)-10s) %(message)s')

### Constants
class fixed:
    D = 3               # dimension of velocity space
    T_B = 1             # temperature of the plates
    L = 0.5             # width of computational domain
    sqr_a = .5          # sound velocity in LBM

args.width = eval(args.width)*args.kn
zeros, ones = np.zeros(fixed.D), np.ones(fixed.D)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
hodge[0, 2, 1] = hodge[2, 1, 0] = hodge[1, 0, 2] = 1
e_x, e_y, e_z = np.eye(fixed.D)

### Auxiliary functions (d = domain, m = model)
delta_y0 = lambda d: d.L/d.N if d.q == 1 else d.L*(1-d.q)/(1-d.q**d.N)
delta_y = lambda d: delta_y0(d) * d.q**np.arange(d.N)
cells = lambda d: d.y0 + np.cumsum(delta_y(d) + np.roll(np.append(delta_y(d)[:-1], 0), 1))/2
half_M = lambda m: -np.sum(m.xi()[...,1] * m.Maxw(Macro(vel=args.U*e_x/2, temp=fixed.T_B)) * (m.xi()[...,1] < 0))
_to_arr = lambda vel: np.atleast_2d(vel)
to_arr = lambda macro: macro if hasattr(macro.rho, '__len__') else Macro(*[ np.array([m]) for m in macro._asdict().values() ])
_from_arr = lambda vel: slice(None) if len(vel.shape) > 1 else 0
from_arr = lambda macro: slice(None) if hasattr(macro.rho, '__len__') else 0

def poly_correction(moments, energy_moments, delta, xi):
    m0, m1_i, m2, m2_ij, m3_i, m4 = moments()
    M2, M3_i, M4 = energy_moments()
    _ = lambda v: v.reshape((-1,1))
    __ = lambda v: v.reshape((-1,fixed.D,1))
    A = np.hstack((
        np.hstack(( _(m0), m1_i, _(m2) )).reshape((-1,1,2+fixed.D)),
        np.dstack(( __(m1_i), m2_ij, __(m3_i) )),
        np.hstack(( _(M2), M3_i, _(M4) )).reshape((-1,1,2+fixed.D))
    ))
    alpha = np.linalg.solve(A, delta)
    sqr_xi = np.einsum('il,il->i', xi, xi)
    psi = np.hstack(( _(np.ones_like(sqr_xi)), xi, _(sqr_xi) ))
    return np.einsum('al,il->ai', alpha, psi)

def poly_correct_moments(model, F, mask, M0, M1_i, M2):
    m0, m1_i, m2, m2_ij, m3_i, m4 = calc_first_moments0(model, F, mask, all=True)
    _ = lambda v: v.reshape((-1,1))
    return poly_correction(
        lambda: (m0, m1_i, m2, m2_ij, m3_i, m4),
        lambda: (m2, m3_i, m4),
        np.hstack(( _(M0), M1_i, _(M2) )),
        model.xi()[mask,:]
    )[_from_arr(F)]

def poly_correct(f, macro, xi):
    # see Aristov, Tcheremissine 1980
    m0, m1_i, m2, m2_ij, m3_i, m4 = calc_first_moments0(None, f, all=True, xi=xi)
    sqr_xi = np.einsum('il,il->i', xi, xi)
    c = np.einsum('il,a->ail', xi, np.ones_like(macro.rho)) - np.einsum('i,al->ail', np.ones_like(sqr_xi), macro.vel)
    sqr_c = np.einsum('ail,ail->ai', c, c)
    M2 = np.einsum('ai,ai->a', f, sqr_c)
    M3_i = np.einsum('ai,ai,il->al', f, sqr_c, xi)
    M4 = np.einsum('ai,ai,i->a', f, sqr_c, sqr_xi)
    _ = lambda v: v.reshape((-1,1))
    exact = np.einsum('a,al->al', macro.rho, np.hstack(( _(np.ones_like(macro.rho)), macro.vel, fixed.D*_(macro.temp)/2. )))
    real = np.hstack(( _(m0), m1_i, _(M2) ))
    return f * (1 + poly_correction(
        lambda: (m0, m1_i, m2, m2_ij, m3_i, m4),
        lambda: (M2, M3_i, M4),
        exact - real, xi
    ))

def entropy_correct(f, macro, xi):
    # see Mieussens 2000
    raise NotImplementedError()

def lagrange_correct(f, macro, xi):
    # see Gamba, Tharkabhushanam 2009
    raise NotImplementedError()

correct_vdf = lambda func, macro, xi: {
    'none': lambda f, m, xi0: f,
    'poly': poly_correct,
    'entropy': entropy_correct,
    'lagrange': lagrange_correct
}[args.correction](func(to_arr(macro)), to_arr(macro), xi)[from_arr(macro)]

vgrid_params = {
    'polynomial': { 'w_min': args.w_min, 'p': args.poly },
    'geometric': { 'q': args.ratio },
}

### DVM velocity grid
def dvm_grid():
    # 1d packing, symm=1 for antisymmetry, symm=-1 for specular reflection
    def dvm_ball(symm):
        _xi_y = _xi[...,1]
        p = np.where((_xi_y > 0)*(_sqr_xi <= args.radius**2))
        n = np.where((_xi_y < 0)*(_sqr_xi <= args.radius**2))
        return np.hstack((n[0], p[0][::symm])), np.hstack((n[1], p[1])), np.hstack((n[2], p[2]))

    ygrid = getattr(vgrid, args.grid.capitalize())(args.radius, args.M*args.My, **vgrid_params.get(args.grid, {}))
    grid = getattr(vgrid, 'Hermite')(args.radius, args.M)
    _Xi = lambda a: ygrid.x if a==1 else grid.x
    _W = lambda a: ygrid.w if a==1 else grid.w
    _I = lambda a: np.ones(_W(a).size)
    _G = lambda a: np.roll((_Xi(a), _I(a), _I(a)), a, axis=0)
    _xi = np.einsum('li,lj,lk->ijkl', _G(0), _G(1), _G(2))
    _sqr_xi = np.einsum('ijkl,ijkl->ijk', _xi, _xi)
    _ball = dvm_ball(1)
    _w = np.einsum('i,j,k->ijk', _W(0), _W(1), _W(2))[_ball]

    I = lambda v, a: np.ones((v.shape[0], _W(a).size))
    G = lambda v, a: np.roll((np.einsum('a,i', np.ones(v.shape[0]), _Xi(a))
        - np.einsum('a,i', v[:,a], _I(a)), I(v, a), I(v, a)), a, axis=0)

    xi = lambda v: np.einsum('lai,laj,lak->aijkl', G(v, 0), G(v, 1), G(v, 2))[(slice(None),) + _ball]
    sqr_xi = lambda v: np.einsum('ail,ail->ai', xi(v), xi(v))
    Maxw = lambda m: np.einsum('a,i,ai->ai', m.rho*(np.pi*m.temp)**(-fixed.D/2), _w,
            np.exp(np.einsum('ai,a->ai', -sqr_xi(m.vel), 1./m.temp)))

    G_tau = lambda m: np.einsum('ail,aim,an,lmn,a->ai', xi(m.vel), xi(m.vel), m.tau, hodge, 1/m.rho/m.temp**2)
    G_qflow = lambda m: 4 * np.einsum('ail,al,ai,a->ai', xi(m.vel), m.qflow,
            np.einsum('ai,a->ai', sqr_xi(m.vel), 1/m.temp/(fixed.D+2)) - .5, 1/m.rho/m.temp**2)
    Grad13 = lambda m: Maxw(m) * (1 + G_tau(m) + G_qflow(m))

    return Model(
        info = 'DVM: ({}x{}x{})'.format(*_sqr_xi.shape),
        weights = _w,
        xi = lambda vel=zeros: xi(_to_arr(vel))[_from_arr(vel)],
        Maxw = lambda macro: correct_vdf(Maxw, macro, _xi[_ball]),
        Grad13 = lambda macro: correct_vdf(Grad13, macro, _xi[_ball])
    )

### LBM velocity grid
def lbm_d3(order, stretch, nodes, weights, full=True):
    lattice = []
    def sort(item):
        return 10*np.sum(item[0]) + np.sum(item[0]/10**np.arange(fixed.D))
    def add_symm(i, node, lattice):
        if i < fixed.D:
            s = np.ones(fixed.D); s[i] *= -1
            if node[i]:
                add_symm(i+1, node, lattice)
                add_symm(i+1, s*node, lattice)
            else:
                add_symm(i+1, node, lattice)
        else:
            lattice += [ (node, weights[0]) ]
    def add_cycle(node, lattice):
        _node = np.array(node)
        while True:
            add_symm(0, _node*stretch*np.sqrt(fixed.sqr_a), lattice)
            _node = np.roll(_node, 1)
            if (_node == node).all():
                break
    for node in nodes:
        add_cycle(node, lattice)
        if len(set(node)) == fixed.D and full:
            add_cycle(node[::-1], lattice)
        weights = np.roll(weights, -1)
    xi, w = np.transpose(sorted(lattice, key=sort))
    return Lattice(np.vstack(xi), np.hstack(w), order)

def lbm_grid():
    _xi, _w, order = lattices[args.lattice].xi, lattices[args.lattice].w, lattices[args.lattice].order
    a, D = fixed.sqr_a, fixed.D
    xi_v = lambda v: np.einsum('il,al', _xi, v)/a
    sqr_xi = np.einsum('il,il->i', _xi, _xi)/a
    sqr = lambda v: np.einsum('al,al,i->ai', v, v, np.ones_like(_w))/a
    weighted = lambda f: np.einsum('i,ai->ai', _w, f)
    weighted_rho = lambda m, f: np.einsum('a,i,ai->ai', m.rho, _w, f)
    xi = lambda v: np.einsum('il,a', _xi, np.ones(v.shape[0])) - np.einsum('i,al', np.ones_like(_w), v)
    Tn = lambda n: 1./np.math.factorial(n) * np.heaviside(order-2*n, 1)

    T1 = lambda m: xi_v(m.vel) * Tn(1)
    T2 = lambda m: ( (xi_v(m.vel))**2 - sqr(m.vel) + np.einsum('a,i', m.temp-1, sqr_xi-D) ) * Tn(2)
    T3 = lambda m: xi_v(m.vel)*( xi_v(m.vel)**2 - 3*sqr(m.vel) + 3*np.einsum('a,i', m.temp-1, sqr_xi-D-2) ) * Tn(3)
    Maxw = lambda m: weighted_rho(m, 1 + T1(m) + T2(m) + T3(m))
    G2 = lambda m: np.einsum('il,im,an,lmn->ai', _xi, _xi, m.tau, hodge)/a * Tn(2)
    G3 = lambda m: xi_v(m.qflow) * (sqr_xi/(fixed.D+2)-1) + G2(m)*xi_v(m.vel) \
        - np.einsum('il,am,an,lmn->ai', _xi, m.vel, m.tau, hodge)/a
    Grad13 = lambda m: Maxw(m) + weighted(G2(m) + G3(m))
    adapter = lambda func, macro: func(to_arr(macro))[from_arr(macro)]

    return Model(
        info = 'LBM: %s' % args.lattice,
        weights = _w,
        xi = lambda vel=zeros: xi(_to_arr(vel))[_from_arr(vel)],
        Maxw = lambda macro: adapter(Maxw, macro),
        Grad13 = lambda macro: adapter(Grad13, macro)
    )

def splot(model, f):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    xi = model.xi()
    m = ( xi[:,2] <= (model.weights0[0])**(1./3)/2 ) * ( xi[:,2] >= 0 )
    #ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.weights0())[m])
    ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.Maxw(Macro()))[m])
    plt.show()

def plot_profiles(solution):
    plt.clf()
    y, h, m = calc_macro(solution)
    U, factor = args.U, args.qflow_factor
    plt.title(domains[0].model.info, loc='left')
    plt.title(domains[1].model.info, loc='right')
    plt.title('k = %g, U = %g' % (args.kn, args.U))
    plt.axvline(x=domains[1].y0, color='k')
    if args.plot_norms:
        for name in norms.keys():
            y, norm = calc_norm(solution, name)
            plt.plot(y, norm, '*-', label=name)
        plt.plot(y, -m.qflow[:,0], 'b*-', label='heat flow')
        plt.semilogy()
        plt.legend()
    else:
        Y, Vel, Tau = np.loadtxt('k0.1.txt').T
        plt.plot(Y, Vel, 'rD--', Y, -2*Tau, 'gD--')
        Y, Vel, Tau, Qflow = np.loadtxt('k0.1-my.txt').T
        plt.plot(Y, factor*Qflow, 'b--')
        plt.plot(y, m.vel[:,0]/U, 'r*-', label='velocity/U')
        plt.plot(y, -m.tau[:,2]/U, 'g*-', label='share stress/U')
        plt.plot(y, -factor*m.qflow[:,0]/U, 'b*-', label='%g*qflow_x/U' % factor)
        plt.plot(y, (m.temp-1)/U**2, 'y', label='(temp-1)/U^2')
        plt.plot(y, m.qflow[:,1]/U**2/args.kn, 'c', label='qflow_y/U^2/k')
        y, norm = calc_norm(solution, args.norm)
        factor = -int(factor*m.qflow[-1,0]/U/norm[-1])
        plt.plot(y, norm*factor, 'k*-', label='%g*%s' % (factor, args.norm))
        if args.plot_ns:
            plt.plot(Y, Y, 'r:')
            plt.plot(Y, args.kn+Y*0, 'g:')
            plt.plot(Y, .1-.4*Y**2, 'y:')
            plt.plot(Y, .8*Y, 'c:')
        legend = plt.legend(loc='upper center')
    plt.xlim(0, .5)
    plt.grid()
    plt.show()
    plt.pause(1e-3)

def calc_macro0(model, F):
    f = _to_arr(F)
    #TODO: optimize using generators
    rho = np.einsum('ai->a', f)
    vel = np.einsum('ai,il,a->al', f, model.xi(), 1./rho)
    c = model.xi(vel)
    sqr_c = np.einsum('ail,ail->ai', c, c)
    csqr_c = np.einsum('ai,ail->ail', sqr_c, c)
    cc_ij = np.einsum('ail,aim->ailm', c, c)
    cc = np.einsum('ailm,lmn', cc_ij, hodge)
    temp = np.einsum('ai,ai,a->a', f, sqr_c, 2./rho)/fixed.D
    qflow = np.einsum('ai,ail->al', f, csqr_c)
    tau = np.einsum('ai,ail->al', f, cc)
    return Macro(*[ m[_from_arr(F)] for m in (rho, vel, temp, tau, qflow) ])

def calc_moment0(model, F, n):
    f = _to_arr(F)
    xi = np.sqrt(np.einsum('il,il->i', model.xi(), model.xi()))
    return np.einsum('ai,i->a', f, xi**n)

def calc_first_moments0(model, F, mask=slice(None), all=False, xi=None):
    f = _to_arr(F)[:,mask]
    if xi is None:
        xi = model.xi()[mask,:]
    sqr_xi = np.einsum('il,il->i', xi, xi)
    m0 = np.einsum('ai->a', f)
    m1_i = np.einsum('ai,il->al', f, xi)
    m2 = np.einsum('ai,i->a', f, sqr_xi)
    if not all:
        return m0, m1_i, m2
    m2_ij = np.einsum('ai,il,im->alm', f, xi, xi)
    m3_i = np.einsum('ai,i,il->al', f, sqr_xi, xi)
    m4 = np.einsum('ai,i,i->a', f, sqr_xi, sqr_xi)
    return m0, m1_i, m2, m2_ij, m3_i, m4

def hermite_polynomials(model):
    xi, I = model.xi() / np.sqrt(fixed.sqr_a), np.eye(fixed.D)
    S = lambda n: 2**(n/2)
    h0 = np.ones_like(xi[:,0]) * S(0)
    h1 = xi * S(1)
    h2 = (np.einsum('il,im->ilm', xi, xi) - np.einsum('i,lm', h0, I)) * S(2)
    h3 = (np.einsum('il,im,in->ilmn', xi, xi, xi) - np.einsum('il,mn', xi, I)
        - np.einsum('im,ln', xi, I) - np.einsum('in,lm', xi, I)) * S(3)
    return h0, h1, h2, h3

def hermite_reconstruct(old_model, new_model, f):
    macro = calc_macro0(old_model, f)
    f = _to_arr(f)
    H_old = hermite_polynomials(old_model)
    H_new = hermite_polynomials(new_model)
    Idx = [ '', 'l', 'lm', 'lmn' ]
    C = lambda n: fixed.sqr_a**n/np.math.factorial(n)
    result = np.zeros((f.shape[0], H_new[0].size))
    for n, (h_old, h_new, idx) in enumerate(zip(H_old, H_new, Idx)):
        a = np.einsum('ai,i{0}'.format(idx), f, h_old)
        result += np.einsum('a{0},i{0}'.format(idx), a, h_new)*C(n)
    return np.einsum('i,ai->ai', new_model.Maxw(Macro()), result)[from_arr(macro)]

def reconstruct(old_model, new_model, f):
    return {
        'grad13': lambda: new_model.Grad13(calc_macro0(old_model, f)),
        'hermite40': lambda: hermite_reconstruct(old_model, new_model, f),
    }[args.reconstruction]()

def calc_macro(solution):
    y, h, rho, vel, temp, tau, qflow = [], [], [], [], [], [], []
    for d, s in zip(domains, solution):
        y += [ cells(d) ]
        h += [ delta_y(d) ]
        for m, m0 in zip((rho, vel, temp, tau, qflow), calc_macro0(d.model, s.f)._asdict().values()):
            m.append(m0)
    s = lambda m: np.hstack(m)
    v = lambda m: np.vstack(m)
    return s(y), s(h), Macro(s(rho), v(vel), s(temp), v(tau), v(qflow))

def calc_norm(solution, norm):
    y, result = [], []
    for d, s in zip(domains, solution):
        y += [ cells(d) ]
        macro = calc_macro0(d.model, s.f)
        maxw, grad13 = d.model.Maxw(macro), d.model.Grad13(macro)
        result += [ norms[norm](s.f, grad13)/norms[norm](s.f, 0) ]
    s = lambda m: np.hstack(m)
    return s(y), s(result)

def calc_moment(solution, n):
    y, result = [], []
    for d, s in zip(domains, solution):
        y += [ cells(d) ]
        result += [ calc_moment0(d.model, s.f, n) ]
    s = lambda m: np.hstack(m)
    return s(y), s(result)

def total_values(domains):
    y, h, macro = calc_macro(domains)
    print('Total mass =', 2*sum(h*macro.rho))
    print('Total half momentum =', np.einsum('a,al', h*macro.rho, macro.vel))
    print('Total energy =', 2*sum(h*macro.rho*(fixed.D*macro.temp/2 + np.einsum('al,al->a', macro.vel, macro.vel))))

def check(f, flux=False):
    if np.sum(np.isnan(f)) > 0:
        logging.error('NaN:')
        print(f)
        raise NameError("NaN has been found!")
    if np.sum(f<0) > 0 and not flux:
        logging.error('Negative:')
        print(np.argwhere(f<0))
        print(f[f<0])
        raise NameError("Negative value has been found!")

# Second-order TVD scheme
def transport(domain, bc, solution, delta_t):
    np.seterr(all='raise') # for all transport threads
    logging.debug('Starting transfer')
    N, xi_y, _xi_y = domain.N, domain.model.xi(np.zeros((domain.N, fixed.D)))[...,1], domain.model.xi()[...,1]
    gamma, h = _xi_y * delta_t,  delta_y(domain)
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    def calc_F(h, f, idx, F, Idx, mask, low_order=False):
        logging.debug(' - calc_F')
        g = np.einsum('i,a', np.abs(gamma), 1/h[idx])
        d1, d2 = f[idx+1] - f[idx], f[idx] - f[idx-1]
        h1, h2 = (h[idx+1] + h[idx])/2, (h[idx] + h[idx-1])/2
        D = lambda d, h: np.einsum('ai,a->ai', np.abs(d), 1/h)[mask]
        DD = lambda d, h, C: np.einsum('ai,a,ai->ai', np.abs(d), 1/h, C)[mask]
        H = np.einsum('a,ai->ai', h[idx], np.sign(d1))[mask]
        small = np.finfo(float).eps
        lim = {
            'none': lambda: 0,
            'lax-wendroff': lambda: D(d1, h1),
            'beam-warming': lambda: D(d2, h2),
            'fromm': lambda: D(d1+d2, h1+h2),
            'minmod': lambda: np.minimum(D(d1, h1), D(d2, h2)),
            'mc': lambda: np.minimum(D(d1+d2, h1+h2), 2*np.minimum(D(d1, h1), D(d2, h2))),
            'wide-third': lambda: np.minimum((DD(d1, h1, 2-g) + DD(d2, h2, 1+g))/3,
                np.minimum(DD(d1, h1, 2/(1-g)), DD(d2, h2, 2/(g + small)))),
            'superbee': lambda: np.maximum(np.minimum(2*D(d1, h1), D(d2, h2)), np.minimum(D(d1, h1), 2*D(d2, h2)))
        }[args.limiter]() * (0 if low_order else 1)
        F[Idx][mask] = f[idx][mask] + (1-g[mask])*H/2 * np.where(d1[mask]*d2[mask] > 0, lim, 0)
    def calc2N(h, f, F, sgn, bc):
        logging.debug(' - calc2N %+d', sgn)
        m0, m1, m2, mN = _mask(sgn), mask(sgn, 1), mask(sgn, 2), mask(sgn, N-2)
        f3[:2], h3[:2] = f[-2:], h[-2:]
        h3[2] = bc.last_cell(h3, f3, m0)
        calc_F(h3, f3, np.array([1]), F, slice(-1, None), m1)   # last flux
        bc.last_flux_is_calculated()
        calc_F(h, f, np.arange(1, N-1), F, slice(2, -1), mN)    # interior fluxes
        check(F[2:][mask(sgn, N-1)], flux=True) # after calc2N
    def calc01(h, f, F, sgn, bc):
        logging.debug(' - calc01 %+d', sgn)
        m0, m1, m2 = _mask(sgn), mask(sgn, 1), mask(sgn, 2)
        bc.first_flux(h[0], f[0], F, m0, calc_F)                # first flux
        f3[1:], h3[1:] = f[:2], h[:2]
        h3[0] = bc.first_cell(h3, f3, F[0], np.abs(gamma)/h[0], m0)
        calc_F(h3, f3, np.array([1]), F, slice(1, 2), m1)       # second flux
        check(F[:2][m2], flux=True) # after calc01
    f, F, f3, h3 = solution.f, solution.F, solution.f3, np.empty(3)
    #F.fill(np.NaN); f3.fill(np.NaN);                            # for debug
    calc2N(h, f, F, 1, bc[1])
    calc2N(h[::-1], f[::-1], F[::-1], -1, bc[0])
    calc01(h, f, F, 1, bc[0])
    calc01(h[::-1], f[::-1], F[::-1], -1, bc[1])
    check(f3, flux=True)
    check(F, flux=True)
    for b in bc:
        b.all_fluxes_are_calculated()
    for b in bc:
        b.wait_for_update()
    for n, b in enumerate(bc):
        b.correct_first_flux(_xi_y, F[-n], _mask(1-2*n))
    f -= np.einsum('i,a,ai->ai', gamma, 1/h, F[1:] - F[:-1])
    #print('bflux for y0=%.2f: %.8e, %.8e' % (domain.y0, np.einsum('i,i' , gamma, F[0]), np.einsum('i,i' , gamma, F[-1])))
    check(f) # after transport

def bgk(domain, bc, solution, delta_t):
    logging.debug('Starting BGK')
    macro = calc_macro0(domain.model, solution.f)
    nu = lambda rho: rho / args.kn
    M = domain.model.Maxw(macro)
    solution.f[:] = np.einsum('ai,a->ai', solution.f - M, np.exp(-delta_t*nu(macro.rho))) + M
    check(solution.f) # after BGK

# Boundary conditions
class Boundary(object):
    def __init__(self, n):
        pass
    def __call__(self, F0, mask):                   # boundary condition, F0 = f{1/2}
        pass
    def _new_f3(self, h, f):
        return np.tile(f, (3,1)), np.tile(h, 3)
    def _calc_flux(self, F, mask, calc_F, h3, f3, low_order=False):
        calc_F(h3, f3, np.array([1]), F, slice(1), np.array([mask]), low_order)
    def last_flux_is_calculated(self):
        pass
    def correct_first_flux(self, xi_y, F, mask):
        pass
    def all_fluxes_are_calculated(self):
        pass
    def wait_for_update(self):
        pass
    def connect(self, solution):
        pass
    def last_cell(self, h3, f3, mask):              # last ghost cell, f3 = f{M-1, M, ghost}
        f3[2][mask] = (2*f3[1] - f3[0])[mask]       # extrapolation by default
        return h3[0]
    def first_flux(self, h1, f1, F, mask, calc_F):  # inflow flux, f1 = f{1}
        F[0][mask] = self(F[0], mask)
    def first_cell(self, h3, f3, F0, gamma1, mask): # first ghost cell, f3 = f{ghost, 1, 2}
        f3[0][mask] = self(f3[1], mask)
        return h3[1]

class Symmetry(Boundary):
    def __call__(self, F0, mask):
        return F0[::-1][mask]
    def last_cell(self, h3, f3, mask):
        f3[2][mask] = self(f3[1], mask)
        return h3[1]

class Diffuse(Boundary):
    def __init__(self, n):
        self._model = domains[n].model
        self._xi_y = self._model.xi()[...,1]
        self._half_M = half_M(self._model)
    def __call__(self, F0, mask):
        rho = np.sum(self._xi_y[::-1][mask] * F0[::-1][mask]) / self._half_M
        return self._model.Maxw(Macro(rho, args.U*e_x/2, fixed.T_B))[mask]
    def first_cell(self, h3, f3, F0, gamma1, mask):
        f3[0][mask] = ((2*F0 - (1-gamma1)*f3[1])/(1+gamma1))[mask]
        return h3[1]

class Couple(Boundary):
    def __init__(self, n, n_partner, idx_partner):
        self._n = n
        self._n_partner = n_partner
        self._idx_partner = idx_partner
        self._lock_F = threading.Event()
        self._lock_f = threading.Event()
    def _partner_cell(self, f, mask):
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        f_partner = self._f_partner[self._idx_partner]
        if model == model_partner:
            f[mask] = f_partner[mask]
        else:
            f[mask] = reconstruct(model_partner, model, f_partner)[mask]
        return delta_y(domains[self._n_partner])[self._idx_partner]
    def connect(self, solution):
        self._f_partner = solution[self._n_partner].f
        self._F_partner = solution[self._n_partner].F
    def first_flux(self, h1, f1, F, mask, calc_F):
        self._lock_F.wait()
        self._lock_F.clear()
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        if model == model_partner:
            # take the prepared flux from the partner
            F[0][mask] = self._F_partner[self._idx_partner][mask]
        else:
            # reconstruct a flux from the partner distribution function (2 cells)
            f3, h3 = self._new_f3(h1, f1)
            idx = 1 + 2*self._idx_partner # -1 or 1
            f3[:2] = reconstruct(model_partner, model, self._f_partner[::-idx][-2:])
            h3[:2] = delta_y(domains[self._n_partner])[::-idx][-2:]
            self._calc_flux(F, mask, calc_F, h3, f3)
    def first_cell(self, h3, f3, F0, gamma1, mask):
        return self._partner_cell(f3[0], mask)
    def last_cell(self, h3, f3, mask):
        return self._partner_cell(f3[2], mask)
    def correct_first_flux(self, xi_y, F, mask):
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        if model != model_partner:
            # conservative correction of the reconstructed flux
            if self._idx_partner:               # correct for both domains
                xi_y_partner = model_partner.xi()[...,1]
                m0, m1_i, m2 = calc_first_moments0(model_partner, xi_y_partner*self._F_partner[self._idx_partner])
                M0, M1_i, M2 = calc_first_moments0(model, xi_y*F)
                logging.debug('Couple BC velocity correction: mass %g, momentum %g, energy %g'
                    % ((m0-M0)[0], (m1_i-M1_i)[0,0], (m2-M2)[0]))
                F[mask] *= 1 + poly_correct_moments(model, xi_y*F, mask, m0-M0, m1_i-M1_i, m2-M2)
    def last_flux_is_calculated(self):
        bcs[self._n_partner][self._idx_partner]._lock_F.set()
    def all_fluxes_are_calculated(self):
        bcs[self._n_partner][self._idx_partner]._lock_f.set()
    def wait_for_update(self):
        self._lock_f.wait()
        self._lock_f.clear()

def create_bcs():
    bcs = []
    for n, boundary in enumerate(boundaries):
        # for each domain
        bc = []
        for p in boundary:
            # for each boundary of domain
            constructor, kwargs = p
            kwargs['n'] = n
            bc.append(constructor(**kwargs))
        bcs.append(bc)
    return bcs

def new_solution(initial):
    solution = [ Solution(d, initial) for d in domains ]
    for bc in bcs:
        for b in bc:
            b.connect(solution)
    return solution

def run_threads(solution, operator, delta_t):
    threads = []
    for d, s, b in zip(domains, solution, bcs):
        thread = threading.Thread(target=operator, args=(d, b, s, delta_t),
            name='%s for domain: y0=%.2f, N=%d' % (operator.__name__, d.y0, d.N))
        thread.start()
        threads.append(thread)
    for t in threads:
        t.join()            # wait until all threads terminate

def solve_bgk():
    U, U0 = args.U, .9*args.U if args.kn < 1 else 0
    solution = new_solution(lambda y: (1+0*y, U0*np.einsum('a,l', y, e_x), 1+0*y))
    y, h, m0 = calc_macro(solution)
    delta_t = min(h)/args.radius/args.step
    total_values(solution)
    if args.plot:
        plt.ion()
        plot_profiles(solution)
    for i in range(args.end):
        # Symmetry second-order splitting
        run_threads(solution, transport, delta_t)
        run_threads(solution, bgk, 2*delta_t)
        run_threads(solution, transport, delta_t)
        [ check(s.f) for s in solution ]    # to stop iteration process in case of negative values
        if args.plot and not i % args.plot:
            y, h, m = calc_macro(solution)
            rho_disp, pxy_mean = sum(h*(m.rho-m0.rho)**2), np.sum(h*m.tau[:,2])
            pxy_disp = sum(h*(m.tau[:,2]-pxy_mean)**2)
            print('%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel_x[-1]/U = %g, qflow_x[-1]/U = %g'
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, m.vel[-1,0]/U, m.qflow[-1,0]/U ))
            plot_profiles(solution)
    total_values(solution)
    y, h, m = calc_macro(solution)
    if args.plot_norms:
        names = list(norms.keys())
        result = [ calc_norm(solution, name)[1] for name in norms.keys() ]
    else:
        names = [ 'vel_x', 'p_xy', 'q_x' ]
        result = [ m.vel[:,0], m.tau[:,2], m.qflow[:,0] ]
    np.savetxt(sys.stdout, np.transpose([y] + result), fmt='%1.5e',
        header='%11s'*(len(names) + 1) % tuple(['y'] + names))
    #splot(domains[-1].model, solution[-1].f[-1])
    #splot(domains[0].model, solution[0].f[0])
    if args.plot:
        plt.ioff(); plt.show()

def tests():
    fmt_float, fmt_str = '%+.2e', '%9s'
    dev2str = lambda x: fmt_str % '0' if np.fabs(x) < 1e1*np.finfo(float).eps else fmt_float % x
    def is_almost_equal(a, b):
        if args.epsilon:
            assert np.sum(np.abs(a-b)) < args.epsilon
    def print_sample(m, m0):
        s = lambda name: m._asdict()[name] - m0._asdict()[name]
        v = lambda name, idx: m._asdict()[name][idx] - m0._asdict()[name][idx]
        columns = ( 'rho', 'temp', 'vel_x', 'vel_y', 'p_xy', 'q_x', 'q_y' )
        values = ( s('rho'), s('temp'), v('vel', 0), v('vel', 1), v('tau', 2), v('qflow', 0), v('qflow', 1) )
        print('Variable: ', ''.join([fmt_str + ' ' for i in range(len(columns))]) % columns)
        print('Deviation:', ''.join([fmt_str + ' ' for i in range(len(values))]) % tuple(map(dev2str, values)))
    def check_macro(model, f, expected):
        computed = calc_macro0(model, f)
        print_sample(computed, expected)
        for m, m0 in zip(computed._asdict().values(), expected._asdict().values()):
            is_almost_equal(m, m0)
    o = lambda y: np.ones_like(y)
    delta = args.U/2
    rho, vel, temp, tau, qflow = 1 + delta, delta*e_x, 1 + delta**2, delta*e_z, delta*e_x/2
    print('Test #1: Maxwell distribution (rho=%g, vel_x=%g, temp=%g)' % (rho, vel[0], temp))
    macro = Macro(rho, vel, temp)
    for name, model in models.items():
        print('-- %s model:' % name)
        check_macro(model, model.Maxw(macro), macro)
    print(''.join(['-' for i in range(50)]))
    print('Test #2: Grad distribution (rho=%g, vel_x=%g, temp=%g, p_xy=%g)' % (rho, vel[0], temp, tau[2]))
    macro = Macro(rho, vel, temp, tau)
    for name, model in models.items():
        print('-- %s model:' % name)
        check_macro(model, model.Grad13(macro), macro)
    print(''.join(['-' for i in range(50)]))
    print('Test #3: Grad distribution (rho=%g, vel_x=%g, temp=%g, qflow_x=%g)' % (rho, vel[0], temp, qflow[0]))
    macro = Macro(rho, vel, temp, qflow=qflow)
    for name, model in models.items():
        print('-- %s model:' % name)
        check_macro(model, model.Grad13(macro), macro)

Macro = namedtuple('Macro', 'rho vel temp tau qflow')
Macro.__new__.__defaults__ = (1., zeros, 1., zeros, zeros)
Lattice = namedtuple('Lattice', 'xi w order')
# some constants for non-space filling lattices
q, Q = np.sqrt(15), np.sqrt(5)
r, s, t = np.sqrt((15+q)/2), np.sqrt(6-q), np.sqrt(9+2*q)
a1, b1, c1, d1 = 1.07182071542885, 1.21422495340964, 2.92338002226218, 2.49326392161601
a2, b2, c2, d2, e2, f2 = 1.45239166695403, 1.12267801869052, 1.37910253014295, 2.66422150193135, 2.40894358132745,  4.81788716265490
lattices = {
    ### Classical
        # 3th order
    'd3q8': lbm_d3(3, 1, [(1,1,1)], np.array([.125])),
        # 5th order
    'd3q15': lbm_d3(5, np.sqrt(3), [(0,0,0), (1,0,0), (1,1,1)], np.array([2, 1, .125])/9),
    'd3q19': lbm_d3(5, np.sqrt(3), [(0,0,0), (1,0,0), (1,1,0)], np.array([3, 18, 36])**-1.),
    'd3q27': lbm_d3(5, np.sqrt(3), [(0,0,0), (1,0,0), (1,1,0), (1,1,1)], np.array([8, 2, .5, .125])/27),
    ### Shan, Yuan, Chen 2006
        # almost minimal of 7th order
    'd3q39': lbm_d3(7, np.sqrt(1.5), [(0,0,0), (1,0,0), (1,1,1), (2,0,0), (2,2,0), (3,0,0)], [1./12, 1./12, 1./27, 2./135, 1./432, 1./1620]),
        # minimal of 5th order (icosahedron)
    'd3q13n': lbm_d3(5, 1, [(0,0,0), (np.sqrt((5-Q)/2),np.sqrt((5+Q)/2),0)], [ 2./5, 1./20 ], full=False),
    ### Nie, Shan, Chen 2008
    'd3q121': lbm_d3(9, 1.19697977039307435897239, [(0,0,0), (1,0,0), (1,1,1), (1,2,0), (2,2,0), (3,0,0), (2,3,0), (2,2,2), (1,1,3), (3,3,3)],
        [ 0.03059162202948600642469, 0.09851595103726339186467, 0.02752500532563812386479, 0.00611102336683342432241, 0.00042818359368108406618,
        0.00032474752708807381296, 0.00001431862411548029405, 0.00018102175157637423759, 0.00010683400245939109491, 0.00000069287508963860285 ]),
    ### Feuchter, Schleifenbaum 2016
        # 5th order
    'd3v15': lbm_d3(5, 1.2247448713915890, [(0,0,0), (2,0,0), (1,1,1)], np.array([7, .5, 1])/18),
        # minimal of 7th order
    'd3v38': lbm_d3(7, 0.86602540378443865, [(1,0,0), (2,0,0), (2,2,2), (2,2,0), (6,0,0)],
        [ 6.7724867724867725e-2, 5.5555555555555556e-2, 4.6296296296296296e-3, 1.8518518518518519e-2, 1.7636684303350970e-4 ]),
        # 7th order
    'd3v64': lbm_d3(7, 0.69965342816864754, [(2,0,0), (1,1,1), (5,0,0), (3,3,3), (3,3,0), (3,1,1)],
        [ 5.9646397884737016e-3, 8.0827437008387392e-2, 1.1345266793939999e-3, 9.5680047874015889e-4, 3.9787631334632013e-3, 1.0641080987258957e-2 ]),
    'd3v96': lbm_d3(7, 0.37787639086813054, [(1,1,1), (3,3,3), (3,1,1), (4,4,4), (7,1,1), (6,6,1)],
        [ 1.2655649299880090e-3, 2.0050978770655310e-2, 2.7543347614356814e-2, 4.9712543563172566e-3, 3.6439016726158895e-3, 1.7168180273737716e-3 ]),
    'd3v112': lbm_d3(7, 0.40531852273291520, [(1,1,1), (3,1,1), (4,4,4), (3,2,2), (7,1,1), (5,5,1)],
        [ 3.3503407500643648e-3, 2.8894128958152456e-2, 4.5930345162087793e-3, 4.4163148398082762e-3, 2.3237070220062610e-3, 3.3847240912752922e-3 ]),
        # minimal of 9th order
    'd3v79': lbm_d3(9, 1, [(0,0,0), (1,0,0), (2,0,0), (1,1,1), (2,2,2), (2,2,0), (6,0,0), (3,3,3), (3,1,1)],
        [ 1.0570987654320988e-1, 3.8095238095238095e-2, 2.8645833333333333e-2, 4.9479166666666667e-2, 5.2083333333333333e-4,
        5.2083333333333333e-3, 2.7557319223985891e-6, 9.6450617283950617e-6, 1.3020833333333333e-3 ]),
    ### Stroud 1971 [ correct formulas in Ortiz 2007 (phd thesis), misprint for t in Shan, Yuan, Chen 2006 ]
        # 7th order
    'd3q27n': lbm_d3(7, 1, [(0,0,0), (r,0,0), (s,s,0), (t,t,t)], [ 8*(90+q)/2205, 2*(135-23*q)/15435, (162+41*q)/6174, (783-202*q)/24696]),
    ### Surmas, Ortiz, Philippi 2009
        # almost 9th order (for moment xi^2 xi_a xi_b)
    'd3v33n': lbm_d3(8, 1, [(0,0,0), (a1,0,0), (b1,b1,b1), (c1,0,0), (d1,d1,0)],
        [ 1.69544317872168e-1, 7.53752058968985e-2, 3.90045337112442e-2, 6.86518217744201e-3, 2.08142366598628e-3 ]),
    'd3v59': lbm_d3(8, 1.20288512331026, [(0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,0,0), (2,2,0), (2,2,2), (3,0,0)],
        [ 9.58789162377528e-2, 7.31047082129148e-2, 3.46588971093380e-3, 3.66108082044515e-2,
        1.59235232232060e-2, 2.52480845105094e-3, 7.26968662515159e-5, 7.65879439346840e-4 ]),
        # 9th order
    'd3v53n': lbm_d3(9, 1, [(0,0,0), (a2,0,0), (b2,b2,b2), (c2,d2,0), (e2,e2,e2), (f2,0,0)],
        [ 1.60453547343974e-1, 7.17493485726724e-2, 3.85210759331158e-2, 4.11522633744856e-3, 2.44958972443801e-4, 2.61083127915713e-5 ]),
    'd3v107': lbm_d3(9, 1.07182071542885, [(0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,0,0), (2,1,0), (2,2,0), (3,1,1), (2,2,2), (4,0,0)],
        [ 7.57516860965017e-2, 6.00912802747447e-2, 3.13606906699535e-3, 3.63392812078012e-2, 1.32169332731492e-2,
        4.48492851172950e-3, 2.48755775808342e-3, 6.07432754970149e-4, 4.64179164402822e-4, 4.51928894609872e-5 ]),
}
Model = namedtuple('Model', 'info weights xi Maxw Grad13')
models = {
    'dvm': dvm_grid(),
    'lbm': lbm_grid()
}
Domain = namedtuple('Domain', 'y0 L model N q')
domains = (
    Domain(0, fixed.L-args.width, models[args.model1], args.N1, 1.),
    Domain(fixed.L-args.width, args.width, models[args.model2], args.N2, args.refinement**(-1./(args.N2-1)))
)
boundaries = (
    [ ( Symmetry, {} ),  ( Couple, { 'n_partner': 1, 'idx_partner': 0 }) ],
    [ ( Couple, { 'n_partner': 0, 'idx_partner': -1 }), ( Diffuse, {}) ]
)
bcs = create_bcs()
norms = {
    'L1': lambda f,g: np.einsum('ai->a', np.abs(f-g)),
    'L2': lambda f,g: np.sqrt(np.einsum('ai->a', (f-g)**2)),
    'Linf': lambda f,g: np.max(np.abs(f-g), axis=1)
}

class Solution(object):
    def __init__(self, domain, initial):
        empty = lambda model, N: np.empty((N,) + model.Maxw(Macro()).shape)
        self.f = domain.model.Maxw(Macro(*initial(cells(domain))))  # distribution function
        self.F = empty(domain.model, domain.N + 1)                  # fluxes between cells
        self.f3 = empty(domain.model, 3)                            # ghost + 2 cells

print(''.join(['=' for i in range(50)]))
for m in models.values():
    print('{}, total = {}, min_w = {:.2e}'.format(m.info, m.weights.size, np.min(m.weights)))
print('xi_max = %g, xi_y_type = %s' % (args.radius, args.grid))
print('Kn = %g, U = %g, cells = %d + %d' % (args.kn, args.U, args.N1, args.N2))
print('Model: (antisym)[ %s | %s ](diffuse), limiter = %s' % (args.model1, args.model2, args.limiter))
print('Width:  |<-- %.3f -->|<-- %.3f -->|' % (fixed.L-args.width, args.width))
print('Delta_y:     | %.4f | %.4f -- %.4f |' %
        (delta_y(domains[0])[0], delta_y(domains[1])[0], delta_y(domains[1])[-1]))
print(''.join(['=' for i in range(50)]))

if args.tests:
    tests()
else:
    solve_bgk()

