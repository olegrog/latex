#!/usr/bin/env python
import sys, argparse, threading, logging, traceback
import numpy as np
from functools import partial
from collections import namedtuple
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-U', type=float, default=2e-2, metavar='<float>', help='difference in velocity of plates')
parser.add_argument('-k', '--kn', type=float, default=2e-1/np.pi**.5, metavar='<float>', help='the Knudsen number')
parser.add_argument('-N1', type=int, default=4, metavar='<int>', help='number of cells in the first domain')
parser.add_argument('-N2', type=int, default=4, metavar='<int>', help='number of cells in the second domain')
parser.add_argument('-M', type=int, default=8, metavar='<int>', help='number of points along each axis in the velocity space')
parser.add_argument('-e', '--end', type=int, default=int(1e2), metavar='<int>', help='maximum total time')
parser.add_argument('-w', '--width', type=str, default='1.2*kn', metavar='<expr>', help='width of the second domain')
parser.add_argument('-s', '--step', type=float, default=1.5, metavar='<float>', help='reduce timestep in given times')
parser.add_argument('-R', '--radius', type=float, default=4., metavar='<float>', help='radius of the velocity grid')
parser.add_argument('-r', '--refinement', type=float, default=2., metavar='<float>', help='ratio of maximum and minumum cell width')
parser.add_argument('-m1', '--model1', default='lbm', metavar='dvm|lbm', help='type of velocity model in the first domain')
parser.add_argument('-m2', '--model2', default='dvm', metavar='dvm|lbm', help='type of velocity model in the second domain')
parser.add_argument('-l', '--lattice', default='d3q19', metavar='d3q15|19|27', help='type of velocity lattice')
parser.add_argument('-i', '--limiter', default='mc', metavar='none|minmod|mc|superbee', help='type of limiter')
parser.add_argument('-c', '--correction', default='poly', metavar='none|poly|entropy|lagrange', help='type of correction for construction of discrete Maxwellian')
parser.add_argument('-p', '--plot', type=int, default=10, metavar='<int>', help='plot every <int> steps')
parser.add_argument('-t', '--tests', action='store_true', help='run some tests instead')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

log_level = logging.DEBUG if args.verbose else logging.INFO
logging.basicConfig(level=log_level, format='(%(threadName)-10s) %(message)s')

### Constants
class fixed:
    T_B = 1             # temperature of the plates
    L = 0.5             # width of computational domain
    sqr_a = .5          # sound velocity in LBM
kn = args.kn
args.width = min(eval(args.width), fixed.L/2)
zeros, ones = np.zeros(3), np.ones(3)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
hodge[0, 2, 1] = hodge[2, 1, 0] = hodge[1, 0, 2] = 1
e_x, e_y, e_z = np.eye(3)

### Auxiliary functions (d = domain, m = model)
delta_y0 = lambda d: d.L/d.N if d.q == 1 else d.L*(1-d.q)/(1-d.q**d.N)
delta_y = lambda d: delta_y0(d) * d.q**np.arange(d.N)
cells = lambda d: d.y0 + np.cumsum(delta_y(d) + np.roll(np.append(delta_y(d)[:-1], 0), 1))/2
half_M = lambda m: -np.sum(m.xi()[...,1] * m.Maxw(Macro(vel=args.U*e_x/2, temp=fixed.T_B)) * (m.xi()[...,1] < 0))
_to_arr = lambda vel: np.atleast_2d(vel)
to_arr = lambda macro: macro if hasattr(macro.rho, '__len__') else Macro(*[ np.array([m]) for m in macro._asdict().values() ])
_from_arr = lambda vel: slice(None) if len(vel.shape) > 1 else 0
from_arr = lambda macro: slice(None) if hasattr(macro.rho, '__len__') else 0

### DVM velocity grid
def poly_correct(f, macro, xi):
    # see Aristov, Tcheremissine 1980
    sqr_xi = np.einsum('il,il->i', xi, xi)
    c = np.einsum('il,a->ail', xi, np.ones_like(macro.rho)) - np.einsum('i,al->ail', np.ones_like(sqr_xi), macro.vel)
    sqr_c = np.einsum('ail,ail->ai', c, c)
    m0 = np.einsum('ai->a', f)
    m1_i = np.einsum('ai,il->al', f, xi)
    m2 = np.einsum('ai,i->a', f, sqr_xi)
    m2_ij = np.einsum('ai,il,im->alm', f, xi, xi)
    m3_i = np.einsum('ai,i,il->al', f, sqr_xi, xi)
    M2 = np.einsum('ai,ai->a', f, sqr_c)
    M3_i = np.einsum('ai,ai,il->al', f, sqr_c, xi)
    M4 = np.einsum('ai,ai,i->a', f, sqr_c, sqr_xi)
    _ = lambda v: v.reshape((-1,1))
    __ = lambda v: v.reshape((-1,3,1))
    exact = np.einsum('a,al->al', macro.rho, np.hstack(( _(np.ones_like(macro.rho)), macro.vel, 1.5*_(macro.temp) )))
    real = np.hstack(( _(m0), m1_i, _(M2) ))
    A = np.hstack((
        np.hstack(( _(m0), m1_i, _(m2) )).reshape((-1,1,5)),
        np.dstack(( __(m1_i), m2_ij, __(m3_i) )),
        np.hstack(( _(M2), M3_i, _(M4) )).reshape((-1,1,5))
    ))
    alpha = np.linalg.solve(A, exact - real)
    psi = np.hstack(( _(np.ones_like(sqr_xi)), xi, _(sqr_xi) ))
    return f * (1 + np.einsum('al,il->ai', alpha, psi))

def entropy_correct(f):
    # see Mieussens 2000
    raise NotImplementedError()

def lagrange_correct(f):
    # see Gamba, Tharkabhushanam 2009
    raise NotImplementedError()

def dvm_grid():
    correct = lambda func, macro: {
        'none': lambda f, macro, xi: f,
        'poly': poly_correct,
        'entropy': entropy_correct,
        'lagrange': lagrange_correct
    }[args.correction](func(to_arr(macro)), to_arr(macro), xi0)[from_arr(macro)]

    # 1d packing, symm=1 for antisymmetry, symm=-1 for specular reflection
    def dvm_ball(symm):
        _xi_y = _xi(zeros)[...,1]
        p = np.where((_xi_y > 0)*(_sqr_xi(zeros) <= args.radius**2))
        n = np.where((_xi_y < 0)*(_sqr_xi(zeros) <= args.radius**2))
        return np.hstack((n[0], p[0][::symm])), np.hstack((n[1], p[1])), np.hstack((n[2], p[2]))

    idx, delta_v = lambda M: np.arange(2*M) - M + .5, args.radius / args.M
    _Xi, _I = delta_v * idx(args.M), np.ones(2*args.M)
    dom = lambda v: np.ones(v.shape[0])
    I = lambda v: np.ones((v.shape[0], 2*args.M))
    _G = lambda v, k: np.roll((_Xi-v[k], _I, _I), k, axis=0)
    G = lambda v, k: np.roll((np.einsum('i,a', _Xi, dom(v)) - np.einsum('a,i', v[:,k], _I), I(v), I(v)), k, axis=0)

    _xi = lambda v: np.einsum('li,lj,lk->ijkl', _G(v, 0), _G(v, 1), _G(v, 2))
    _sqr_xi = lambda v: np.einsum('ijkl,ijkl->ijk', _xi(v), _xi(v))
    _ball = dvm_ball(1)
    xi0 = _xi(zeros)[_ball]

    xi = lambda v: np.einsum('lai,laj,lak->aijkl', G(v, 0), G(v, 1), G(v, 2))[(slice(None),) + _ball]
    sqr_xi = lambda v: np.einsum('ail,ail->ai', xi(v), xi(v))
    Maxw = lambda m: np.einsum('a,ai->ai', m.rho/(np.pi*m.temp)**1.5*delta_v**3, np.exp(np.einsum('ai,a->ai', -sqr_xi(m.vel), 1./m.temp)))

    G1 = lambda m: np.einsum('ail,aim,an,lmn,a->ai', xi(m.vel), xi(m.vel), m.tau, hodge, 1/m.rho/m.temp**2)
    G2 = lambda m: .8 * np.einsum('ail,al,ai,a->ai', xi(m.vel), m.qflow, np.einsum('ai,a->ai', sqr_xi(m.vel), 1/m.temp) - 2.5, 1/m.rho/m.temp**2)
    Grad13 = lambda m: Maxw(m) * (1 + G1(m) + G2(m))

    return Model(
        info = 'DVM: (%d)^3' % (2*args.M),
        weights = np.ones_like(_ball[0]) * delta_v**3,
        xi = lambda vel=zeros: xi(_to_arr(vel))[_from_arr(vel)],
        Maxw = lambda macro: correct(Maxw, macro),
        Grad13 = lambda macro: correct(Grad13, macro)
    )

### LBM velocity grid
def lbm_d3(stretch, nodes, weights, full=True):
    lattice = []
    def sort(item):
        return 10*np.sum(item[0]) + np.sum(item[0]/10**np.arange(3.))
    def add_symm(i, node, lattice):
        if i < 3:
            s = np.ones(3); s[i] *= -1
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
        if len(set(node)) == 3 and full:
            add_cycle(node[::-1], lattice)
        weights = np.roll(weights, -1)
    xi, w = np.transpose(sorted(lattice, key=sort))
    return Lattice(np.vstack(xi), np.hstack(w))

def lbm_grid():
    _xi, _w = lattices[args.lattice].xi, lattices[args.lattice].w
    a = fixed.sqr_a
    xi_v = lambda v: np.einsum('il,al', _xi, v)/a
    sqr_xi = np.einsum('il,il->i', _xi, _xi)/a
    sqr = lambda v: np.einsum('al,al,i->ai', v, v, np.ones_like(_w))/a
    weighted = lambda f: np.einsum('i,ai->ai', _w, f)
    weighted_rho = lambda m, f: np.einsum('a,i,ai->ai', m.rho, _w, f)
    xi = lambda v: np.einsum('il,a', _xi, np.ones(v.shape[0])) - np.einsum('i,al', np.ones_like(_w), v)

    T1 = lambda m: xi_v(m.vel)
    T2 = lambda m: ( (xi_v(m.vel))**2 - sqr(m.vel) + np.einsum('a,i', m.temp-1, sqr_xi-3) ) / 2
    T3 = lambda m: xi_v(m.vel)*( xi_v(m.vel)**2 - 3*sqr(m.vel) + 3*np.einsum('a,i', m.temp-1, sqr_xi-5) ) / 6
    Maxw = lambda m: weighted_rho(m, 1 + T1(m) + T2(m) + T3(m))
    G2 = lambda m: np.einsum('il,im,an,lmn->ai', _xi, _xi, m.tau, hodge)
    G3 = lambda m: xi_v(m.qflow) * (sqr_xi/5-1) + G2(m)*xi_v(m.vel) - 2*np.einsum('il,am,an,lmn->ai', _xi, m.vel, m.tau, hodge)
    Grad13 = lambda m: Maxw(m) + weighted(G2(m) + G3(m))
    correct = lambda func, macro: func(to_arr(macro))[from_arr(macro)]

    return Model(
        info = 'LBM: %s' % args.lattice,
        weights = _w,
        xi = lambda vel=zeros: xi(_to_arr(vel))[_from_arr(vel)],
        Maxw = lambda macro: correct(Maxw, macro),
        Grad13 = lambda macro: correct(Grad13, macro)
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
    U, factor = args.U, 40
    Y, Vel, Tau = np.loadtxt('k1e-1.txt').T
    plt.plot(Y, Vel, 'rD--', Y, -2*Tau, 'gD--')                 # k = 0.1
    Y, Vel, Vel2, Tau, Qflow = np.loadtxt('k1e-1-my.txt').T
    plt.plot(Y, factor*Qflow, 'b--')
    #Y, Vel, Tau = np.loadtxt('k1e-0.txt').T
    #plt.plot(Y, Vel, 'rs--', Y, -2*Tau, 'gs--')                 # k = 1
    #plt.plot(Y, Y, 'r-.')                                       # k = 0
    #plt.plot(Y, 0*Y + np.pi**-.5, 'g-.')                        # k = \infty
    plt.plot(y, m.vel[:,0]/U, 'r*-', label='velocity/U')
    plt.plot(y, -m.tau[:,2]/U, 'g*-', label='share stress/U')
    plt.plot(y, -factor*m.qflow[:,0]/U, 'b*-', label='%d*heat flow/U' % factor)
    legend = plt.legend(loc='upper center')
    plt.title(domains[0].model.info, loc='left')
    plt.title(domains[1].model.info, loc='right')
    plt.title('Kn = %g' % args.kn)
    plt.axvline(x=domains[1].y0, color='k')
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
    temp = 2./3*np.einsum('ai,ai,a->a', f, sqr_c, 1./rho)
    qflow = np.einsum('ai,ail->al', f, csqr_c)
    tau = np.einsum('ai,ail->al', f, cc)
    return Macro(*[ m[_from_arr(F)] for m in (rho, vel, temp, tau, qflow) ])

def reconstruct(old_model, new_model, f):
    macro = calc_macro0(old_model, f)
    return new_model.Grad13(macro)

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

def total_values(domains):
    y, h, macro = calc_macro(domains)
    print 'Total mass =', 2*sum(h*macro.rho)

def check(f):
    if np.sum(np.isnan(f)) > 0:
        raise NameError("NaN has been found!")
    if np.sum(f<0) > 0:
        logging.error('Negative:')
        print np.argwhere(f<0)
        print f[f<0]
        raise NameError("Negative value has been found!")

# Second-order TVD scheme
def transport(domain, bc, solution, delta_t):
    np.seterr(all='raise') # for all transport threads
    logging.debug('Starting transfer')
    N, xi_y, _xi_y = domain.N, domain.model.xi(np.zeros((domain.N, 3)))[...,1], domain.model.xi()[...,1]
    gamma, h = _xi_y * delta_t,  delta_y(domain)
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    def calc_F(h, f, idx, F, Idx, mask):
        logging.debug('  - calc_F')
        g = np.einsum('i,a', gamma, 1/h[idx])[mask]
        d1, d2 = f[idx+1] - f[idx], f[idx] - f[idx-1]
        h1, h2 = (h[idx+1] + h[idx])/2, (h[idx] + h[idx-1])/2
        D = lambda d, h: np.einsum('ai,a->ai', np.abs(d), 1/h)[mask]
        H = np.einsum('a,ai->ai', h[idx], np.sign(d1))[mask]
        lim = {
            'none': 0,
            'minmod': np.minimum(D(d1, h1), D(d2, h2)),
            'mc': np.minimum(D(d1+d2, h1+h2), 2*np.minimum(D(d1, h1), D(d2, h2))),
            'superbee': np.maximum(np.minimum(2*D(d1, h1), D(d2, h2)), np.minimum(D(d1, h1), 2*D(d2, h2)))
        }[args.limiter]
        F[Idx][mask] = f[idx][mask] + (1-g)*H/2 * np.where(d1[mask]*d2[mask] > 0, lim, 0)
    def calc2N(h, f, F, sgn, bc):
        logging.debug(' - calc2N %+d', sgn)
        m0, m1, m2, mN = _mask(sgn), mask(sgn, 1), mask(sgn, 2), mask(sgn, N-2)
        f3[:2][m2], h3[:2] = f[-2:][m2], h[-2:]
        h3[2] = bc.last_cell(h3, f3, f, m0)
        calc_F(h3, f3, np.array([1]), F, slice(-1, None), m1)   # last flux
        bc.last_flux_is_calculated()
        calc_F(h, f, np.arange(1, N-1), F, slice(2, -1), mN)    # interior fluxes
        #check(F[2:][mask(sgn, N-1)]) # after calc2N
    def calc01(h, f, F, sgn, bc):
        logging.debug(' - calc01 %+d', sgn)
        m0, m1, m2 = _mask(sgn), mask(sgn, 1), mask(sgn, 2)
        f2, h2 = f[:2], h[:2]
        bc.first_flux(h2, f2, F, m0, calc_F)                    # first flux
        f3[1:][m2], h3[1:] = f2[m2], h2
        h3[0] = bc.first_cell(h3, f3, f, m0)
        calc_F(h3, f3, np.array([1]), F, slice(1, 2), m1)       # second flux
        #check(F[:2][m2]) # after calc01
    f, F, f3, h3 = solution.f, solution.F, solution.f3, np.empty(3)
    #F.fill(np.NaN); f3.fill(np.NaN);                            # for debug
    calc2N(h, f, F, 1, bc[1])
    calc2N(h[::-1], f[::-1], F[::-1], -1, bc[0])
    calc01(h, f, F, 1, bc[0])
    calc01(h[::-1], f[::-1], F[::-1], -1, bc[1])
    #check(f3)
    #check(F)
    for b in bc:
        b.all_fluxes_are_calculated()
    for b in bc:
        b.wait_for_update()
    for n, b in enumerate(bc):
        b.correct_first_flux(_xi_y, F[-n], _mask(1-2*n))
    f -= np.einsum('i,a,ai->ai', gamma, 1/h, F[1:] - F[:-1])
    #print 'bflux for y0=%.2f: %.8e, %.8e' % (domain.y0, np.einsum('i,i' , gamma, F[0]), np.einsum('i,i' , gamma, F[-1]))
    check(f) # after transport

def bgk(domain, bc, solution, delta_t):
    logging.debug('Starting BGK')
    macro = calc_macro0(domain.model, solution.f)
    nu = lambda rho: 2/(np.pi**.5*args.kn) * rho
    M = domain.model.Maxw(macro)
    solution.f[:] = np.einsum('ai,a->ai', solution.f - M, np.exp(-delta_t*nu(macro.rho))) + M
    check(solution.f) # after BGK

# Boundary conditions
class Boundary(object):
    def __init__(self, n):
        pass
    def _new_f3(self, h2, f2):
        return np.tile(f2[0], (3,1)), np.tile(h2[0], 3)
    def _calc_flux(self, F, mask, calc_F, h3, f3):
        calc_F(h3, f3, np.array([1]), F, slice(1), np.array([mask]))
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
    def last_cell(self, h3, f3, f, mask):           # last ghost cell
        f3[2][mask] = (2*f[-1] - f[-2])[mask]       # extrapolation by default
        return h3[0]
    def first_flux(self, h2, f2, F, mask, calc_F):
        self(F[0], F[0], mask)
    def first_cell(self, h3, f3, f, mask):          # first ghost cell
        self(f3[0], f[0], mask)
        return h3[1]

class Symmetry(Boundary):
    def __call__(self, F, F0, mask):
        F[mask] = F0[::-1][mask]
    def first_flux(self, h2, f2, F, mask, calc_F):
        #return self(F[0], F[0], mask)               # uncomment to check mass conservation
        f3, h3 = self._new_f3(h2, f2)
        f3[:2] = f2[::-1,::-1]
        h3[:2] = h2[::-1]
        self._calc_flux(F, mask, calc_F, h3, f3)
    def last_cell(self, h3, f3, f, mask):
        self(f3[2], f[-1], mask)
        return h3[1]
    def correct_first_flux(self, xi_y, F, mask):
        corr = -np.sum((xi_y*F)[~mask]) / np.sum((xi_y*F)[mask])
        logging.debug('Symmetry BC correction %.10e' % corr)
        F[mask] *= corr

class Diffuse(Boundary):
    def __init__(self, n):
        self._model = domains[n].model
        self._xi_y = self._model.xi()[...,1]
        self._half_M = half_M(self._model)
    def __call__(self, F, F0, mask):
        rho = np.sum(self._xi_y[::-1][mask]*F0[::-1][mask])/self._half_M
        F[mask] = self._model.Maxw(Macro(rho, args.U*e_x/2, fixed.T_B))[mask]

class Couple(Boundary):
    def __init__(self, n, n_partner, idx_partner):
        self._n = n
        self._n_partner = n_partner
        self._idx_partner = idx_partner
        self._lock_F = threading.Event()
        self._lock_f = threading.Event()
    def __call__(self, f, mask):
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
    def first_flux(self, h2, f2, F, mask, calc_F):
        self._lock_F.wait()
        self._lock_F.clear()
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        if model == model_partner:
            # take the prepared flux from the partner
            F[0][mask] = self._F_partner[self._idx_partner][mask]
        else:
            # reconstruct a flux from the partner distribution function (2 cells)
            f3, h3 = self._new_f3(h2, f2)
            idx = 1 + 2*self._idx_partner # -1 or 1
            f3[:2] = reconstruct(model_partner, model, self._f_partner[::-idx][-2:])
            h3[:2] = delta_y(domains[self._n_partner])[::-idx][-2:]
            self._calc_flux(F, mask, calc_F, h3, f3)
    def first_cell(self, h3, f3, f, mask):
        return self(f3[0], mask)
    def last_cell(self, h3, f3, f, mask):
        return self(f3[2], mask)
    def correct_first_flux(self, xi_y, F, mask):
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        if model != model_partner:
            # conservative correction of the reconstructed flux
            if self._idx_partner:               # correct for both domains
                xi_y_partner = model_partner.xi()[...,1]
                corr = np.sum(xi_y_partner*self._F_partner[self._idx_partner]) - np.sum(xi_y*F)
                logging.debug('Couple BC velocity correction: %g' % corr)
                #F *= 1 + 2*xi_y*corr           # satisfactory correction for small Ma
                F *= 1 + corr/xi_y/np.sum(F)    # correction for all Ma
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
    plt.ion()
    plot_profiles(solution)
    for i in xrange(args.end):
        # Symmetry second-order splitting
        run_threads(solution, transport, delta_t)
        run_threads(solution, bgk, 2*delta_t)
        run_threads(solution, transport, delta_t)
        [ check(s.f) for s in solution ]
        if not i % args.plot:
            y, h, m = calc_macro(solution)
            rho_disp, pxy_mean = sum(h*(m.rho-m0.rho)**2), np.sum(h*m.tau[:,2])
            pxy_disp = sum(h*(m.tau[:,2]-pxy_mean)**2)
            print '%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel_x[-1]/U = %g, qflow_x[-1]/U = %g' \
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, m.vel[-1,0]/U, m.qflow[-1,0] )
            plot_profiles(solution)
    total_values(solution)
    y, h, m = calc_macro(solution)
    names = ('y', 'vel_x', 'p_xy', 'q_x')
    np.savetxt(sys.stdout, np.transpose((y, m.vel[:,0], m.tau[:,2], m.qflow[:,0])), fmt='%1.4e',
        header='%10s'*len(names) % tuple(names))
    #splot(domains[-1].model, solution[-1].f[-1])
    #splot(domains[0].model, solution[0].f[0])
    plt.ioff(); plt.show()

def tests():
    def is_almost_equal(a, b):
        assert np.sum(np.abs(a-b)) < 1e-5
    def print_sample(m, m0):
        s = lambda name: m._asdict()[name] - m0._asdict()[name]
        v = lambda name, idx: m._asdict()[name][idx] - m0._asdict()[name][idx]
        columns = ( 'rho', 'temp', 'vel_x', 'vel_y', 'p_xy', 'q_x', 'q_y' )
        values = ( s('rho'), s('temp'), v('vel', 0), v('vel', 1), v('tau', 2), v('qflow', 0), v('qflow', 1) )
        print 'Variable: ', ''.join(['%9s ' for i in range(len(columns))]) % columns
        print 'Deviation:', ''.join(['%+.2e ' for i in range(len(values))]) % values
    def check_macro(model, f, expected):
        computed = calc_macro0(model, f)
        print_sample(computed, expected)
        for m, m0 in zip(computed._asdict().values(), expected._asdict().values()):
            is_almost_equal(m, m0)
    o = lambda y: np.ones_like(y)
    delta = args.U/2
    rho, vel, temp, tau, qflow = 1 + delta, delta*e_x, 1 + delta**2, delta*e_z, delta*e_x/2
    print 'Test #1: Maxwell distribution (rho=%g, vel_x=%g, temp=%g)' % (rho, vel[0], temp)
    macro = Macro(rho, vel, temp)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, model.Maxw(macro), macro)
    print ''.join(['-' for i in range(50)])
    print 'Test #2: Grad distribution (rho=%g, vel_x=%g, temp=%g, p_xy=%g)' % (rho, vel[0], temp, tau[2])
    macro = Macro(rho, vel, temp, tau)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, model.Grad13(macro), macro)
    print ''.join(['-' for i in range(50)])
    print 'Test #3: Grad distribution (rho=%g, vel_x=%g, temp=%g, qflow_x=%g)' % (rho, vel[0], temp, qflow[0])
    macro = Macro(rho, vel, temp, qflow=qflow)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, model.Grad13( macro), macro)

Macro = namedtuple('Macro', 'rho vel temp tau qflow')
Macro.__new__.__defaults__ = (1., zeros, 1., zeros, zeros)
Lattice = namedtuple('Lattice', 'xi w')
# some constants for non-space filling lattices
q, Q = np.sqrt(15), np.sqrt(5)
r, s, t = np.sqrt((15+q)/2), np.sqrt(6-q), np.sqrt(9+q)
R, S = np.sqrt((5 + Q)/2), np.sqrt((5-Q)/2)
lattices = {
    'd3q15': lbm_d3(np.sqrt(3), [(0,0,0), (1,0,0), (1,1,1)], np.array([2, 1, .125])/9),
    'd3q19': lbm_d3(np.sqrt(3), [(0,0,0), (1,0,0), (1,1,0)], np.array([3, 18, 36])**-1.),
    'd3q27': lbm_d3(np.sqrt(3), [(0,0,0), (1,0,0), (1,1,0), (1,1,1)], np.array([8, 2, .5, .125])/27),
    # Shan, Yuan, Chen 2006
    'd3q39': lbm_d3(np.sqrt(1.5), [(0,0,0), (1,0,0), (1,1,1), (2,0,0), (2,2,0), (3,0,0)], [1./12, 1./12, 1./27, 2./135, 1./432, 1./1620]),
    # Nie, Shan, Chen 2008
    'd3q121': lbm_d3(1.19697977039307435897239, [(0,0,0), (1,0,0), (1,1,1), (1,2,0), (2,2,0), (3,0,0), (2,3,0), (2,2,2), (1,1,3), (3,3,3)],
        [ 0.03059162202948600642469, 0.09851595103726339186467, 0.02752500532563812386479, 0.00611102336683342432241, 0.00042818359368108406618,
        0.00032474752708807381296, 0.00001431862411548029405, 0.00018102175157637423759, 0.00010683400245939109491, 0.00000069287508963860285 ]),
    # Feuchter, Schleifenbaum 2016
    'd3v15': lbm_d3(1.2247448713915890, [(0,0,0), (2,0,0), (1,1,1)], np.array([7, .5, 1])/18),
    'd3v64': lbm_d3(0.69965342816864754, [(2,0,0), (1,1,1), (5,0,0), (3,3,3), (3,3,0), (3,1,1)],
        [ 5.9646397884737016e-3, 8.0827437008387392e-2, 1.1345266793939999e-3, 9.5680047874015889e-4, 3.9787631334632013e-3, 1.0641080987258957e-2 ]),
    'd3v96': lbm_d3(0.37787639086813054, [(1,1,1), (3,3,3), (3,1,1), (4,4,4), (7,1,1), (6,6,1)],
        [ 1.2655649299880090e-3, 2.0050978770655310e-2, 2.7543347614356814e-2, 4.9712543563172566e-3, 3.6439016726158895e-3, 1.7168180273737716e-3 ]),
    'off-d3q13': lbm_d3(1, [(0,0,0), (R,S,0)], [ 2./5, 1./20 ], full=False),
    # Stroud 1971
    'off-d3q27': lbm_d3(1, [(0,0,0), (r,0,0), (s,s,0), (t,t,t)], [ (720+8*q)/2205, (270-46*q)/15435, (162+41*q)/6174, (783-202*q)/24696])
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

class Solution(object):
    def __init__(self, domain, initial):
        empty = lambda model, N: np.empty((N,) + model.Maxw(Macro()).shape)
        self.f = domain.model.Maxw(Macro(*initial(cells(domain))))  # distribution function
        self.F = empty(domain.model, domain.N + 1)                  # fluxes between cells
        self.f3 = empty(domain.model, 3)                            # ghost + 2 cells

print ''.join(['=' for i in range(50)])
print 'DVM: xi_max = %g, grid=(%d)^3, total = %d' % (args.radius, 2*args.M, models['dvm'].xi().size/3)
print 'LBM: type = %s, total = %d' % (args.lattice, models['lbm'].xi().size/3)
print 'Kn = %g, U = %g, cells = %d + %d' % (args.kn, args.U, args.N1, args.N2)
print 'Model: (antisym)[ %s | %s ](diffuse), limiter = %s' % (args.model1, args.model2, args.limiter)
print 'Width:  |<-- %.3f -->|<-- %.3f -->|' % (fixed.L-args.width, args.width)
print 'Delta_y:     | %.4f | %.4f |' % ((fixed.L-args.width)/args.N1, args.width/args.N2)
print ''.join(['=' for i in range(50)])

if args.tests:
    tests()
else:
    solve_bgk()

