#!/usr/bin/env python
import sys, argparse, threading, logging
import numpy as np
from functools import partial
from collections import namedtuple
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-N1', type=int, default=4, metavar='<int>', help='number of cells in the first domain')
parser.add_argument('-N2', type=int, default=4, metavar='<int>', help='number of cells in the second domain')
parser.add_argument('-M', type=int, default=8, metavar='<int>', help='number of points along each axis in the velocity space')
parser.add_argument('-e', '--end', type=int, default=int(1e2), metavar='<int>', help='maximum total time')
parser.add_argument('-w', '--width', type=str, default='1.2*kn', metavar='<expr>', help='width of the second domain')
parser.add_argument('-s', '--step', type=float, default=1, metavar='<float>', help='reduce timestep in given times')
parser.add_argument('-R', '--radius', type=float, default=4., metavar='<float>', help='radius of the velocity grid')
parser.add_argument('-r', '--refinement', type=float, default=2., metavar='<float>', help='ratio of maximum and minumum cell width')
parser.add_argument('-m1', '--model1', default='lbm', metavar='dvm|lbm', help='type of velocity model in the first domain')
parser.add_argument('-m2', '--model2', default='dvm', metavar='dvm|lbm', help='type of velocity model in the second domain')
parser.add_argument('-l', '--lattice', default='d3q19', metavar='d3q15|19|27', help='type of velocity lattice')
parser.add_argument('-o', '--order', type=int, default=2, metavar='1|2', help='order of approximation')
parser.add_argument('-k', '--kn', type=float, default=2e-1/np.pi**.5, metavar='<float>', help='the Knudsen number')
parser.add_argument('-p', '--plot', type=int, default=10, metavar='<int>', help='plot every <int> steps')
parser.add_argument('-U', type=float, default=2e-2, metavar='<float>', help='difference in velocity of plates')
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
e_x, e_y, e_z = np.eye(3)

### Auxiliary functions (d = domain, m = model)
delta_y0 = lambda d: d.L/d.N if d.q == 1 else d.L*(1-d.q)/(1-d.q**d.N)
delta_y = lambda d: delta_y0(d) * d.q**np.arange(d.N)
cells = lambda d: d.y0 + np.cumsum(delta_y(d) + np.roll(np.append(delta_y(d)[:-1], 0), 1))/2
half_M = lambda m: -np.sum(m.xi0()[...,1] * m.Maxw0(Macro(vel=args.U*e_x/2, temp=fixed.T_B)) * (m.xi0()[...,1] < 0))

### DVM velocity grid
def dvm_grid():
    idx, delta_v = lambda M: np.arange(2*M) - M + .5, args.radius / args.M
    _Xi, _I = delta_v * idx(args.M), np.ones(2*args.M)
    dom = lambda v: np.ones(v.shape[0])
    I = lambda v: np.ones((v.shape[0], 2*args.M))
    _G = lambda v, k: np.roll((_Xi-v[k], _I, _I), k, axis=0)
    G = lambda v, k: np.roll((np.einsum('i,a', _Xi, dom(v)) - np.einsum('a,i', v[:,k], _I), I(v), I(v)), k, axis=0)

    _xi = lambda v: np.einsum('li,lj,lk->ijkl', _G(v, 0), _G(v, 1), _G(v, 2))
    _sqr_xi = lambda v: np.einsum('ijkl,ijkl->ijk', _xi(v), _xi(v))
    _Maxw = lambda m: m.rho/(np.pi*m.temp)**1.5 * np.exp(-_sqr_xi(m.vel)/m.temp)*delta_v**3

    xi = lambda v: np.einsum('lai,laj,lak->aijkl', G(v, 0), G(v, 1), G(v, 2))
    sqr_xi = lambda v: np.einsum('aijkl,aijkl->aijk', xi(v), xi(v))
    Maxw = lambda m: np.einsum('a,aijk->aijk', m.rho/(np.pi*m.temp)**1.5*delta_v**3, np.exp(np.einsum('aijk,a->aijk', -sqr_xi(m.vel), 1./m.temp)))

    # 1d packing, symm=1 for antisymmetry, symm=-1 for specular reflection
    def dvm_ball(symm):
        _xi_y = _xi(zeros)[...,1]
        p = np.where((_xi_y > 0)*(_sqr_xi(zeros) <= args.radius**2))
        n = np.where((_xi_y < 0)*(_sqr_xi(zeros) <= args.radius**2))
        return np.hstack((n[0], p[0][::symm])), np.hstack((n[1], p[1])), np.hstack((n[2], p[2]))
    ball = dvm_ball(1)

    return Model(
        info = 'DVM: (%d)^3' % (2*args.M),
        xi0 = lambda v=zeros: _xi(v)[ball],
        Maxw0 = lambda macro: _Maxw(macro)[ball],
        weights0 = np.ones_like(ball[0]) * delta_v**3,
        xi = lambda v: xi(v)[(slice(None),) + ball],
        Maxw = lambda macro: Maxw(macro)[(slice(None),) + ball]
    )

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
    xi_v = lambda v: np.einsum('il,al', lattices[args.lattice].xi, v)
    sqr_xi = np.einsum('il,il->i', lattices[args.lattice].xi, lattices[args.lattice].xi)
    sqr = lambda v: np.einsum('al,al,i->ai', v, v, lat)
    weighted = lambda m, f: np.einsum('a,i,ai->ai', m.rho, lattices[args.lattice].w, f)
    lat = np.ones_like(lattices[args.lattice].w)
    xi = lambda v: np.einsum('il,a', lattices[args.lattice].xi, np.ones(v.shape[0])) - np.einsum('i,al', lat, v)
    a = fixed.sqr_a
    T1 = lambda m: xi_v(m.vel)/a
    T2 = lambda m: ( (xi_v(m.vel)/a)**2 - sqr(m.vel)/a + np.einsum('a,i', m.temp-1, sqr_xi/a-3) ) / 2
    T3 = lambda m: xi_v(m.vel)/a*( (xi_v(m.vel)/a)**2 - 3*sqr(m.vel)/a + 3*np.einsum('a,i', m.temp-1, sqr_xi/a-5) ) / 6
    Maxw = lambda m: weighted(m, 1 + T1(m) + T2(m) + T3(m))
    return Model(
        info = 'LBM: %s' % args.lattice,
        xi0 = lambda vel=zeros: xi(np.array([vel]))[0],
        Maxw0 = lambda macro: Maxw(Macro(*[ np.array([m]) for m in macro._asdict().values() ]))[0],
        weights0 = lattices[args.lattice].w,
        xi = xi,
        Maxw = Maxw
    )

def splot(model, f):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    xi = model.xi0()
    m = ( xi[:,2] <= (model.weights0[0])**(1./3)/2 ) * ( xi[:,2] >= 0 )
    #ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.weights0())[m])
    ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.Maxw0(Macro()))[m])
    plt.show()

def plot_profiles(solution):
    plt.clf()
    y, h, m = calc_macro(solution)
    U = args.U
    Y, Vel, Tau = np.loadtxt('k1e-1.txt').T
    plt.plot(Y, Vel, 'rD--', Y, -2*Tau, 'gD--')                 # k = 0.1
    Y, Vel, Tau = np.loadtxt('k1e-0.txt').T
    plt.plot(Y, Vel, 'rs--', Y, -2*Tau, 'gs--')                 # k = 1
    plt.plot(Y, Y, 'r-.')                                       # k = 0
    plt.plot(Y, 0*Y + np.pi**-.5, 'g-.')                        # k = \infty
    plt.plot(y, m.vel[:,0]/U, 'r*-', label='velocity/U')
    plt.plot(y, -m.tau[:,2]/U, 'g*-', label='share stress/U')
    plt.plot(y, -m.qflow[:,0]/U/U, 'b*-', label='heat flow/U^2')
    legend = plt.legend(loc='upper center')
    plt.title(domains[0].model.info, loc='left')
    plt.title(domains[1].model.info, loc='right')
    plt.title('Kn = %g' % args.kn)
    plt.axvline(x=domains[1].y0, color='k')
    plt.show()
    plt.pause(1e-3)

def calc_macro0(model, f):
    f = [f] if len(f.shape) == 1 else f
    #TODO: optimize using generators
    rho = np.einsum('ai->a', f)
    vel = np.einsum('ai,il,a->al', f, model.xi0(), 1./rho)
    c = model.xi(vel)
    sqr_c = np.einsum('ail,ail->ai', c, c)
    csqr_c = np.einsum('ai,ail->ail', sqr_c, c)
    cc_ij = np.einsum('ail,aim->ailm', c, c)
    cc = np.einsum('ailm,lmn', cc_ij, hodge)
    temp = 2./3*np.einsum('ai,ai,a->a', f, sqr_c, 1./rho)
    qflow = np.einsum('ai,ail->al', f, csqr_c)
    tau = 2*np.einsum('ai,ail->al', f, cc)
    return Macro(rho, vel, temp, tau, qflow)

def grad13(model, macro):
    if not hasattr(macro.rho, '__len__'):
        macro = Macro(*[ np.array([m]) for m in macro._asdict().values() ])
    c = model.xi(macro.vel)
    sqr_c = np.einsum('ail,ail->ai', c, c)
    #### multiply by 2 due to tau*hodge contains only 3 nonzero elements
    H2 = lambda m: 2 * np.einsum('ail,aim,an,lmn,a->ai', c, c, m.tau, hodge, 1/m.rho/m.temp**2)
    H3 = lambda m: .8 * np.einsum('ail,al,ai,a->ai', c, m.qflow, np.einsum('ai,a->ai', sqr_c, 1/m.temp) - 2.5, 1/m.rho/m.temp**2)
    idx = 0 if len(macro.vel) == 1 else slice(None)
    return (model.Maxw(macro) * (1 + H2(macro) + H3(macro)))[idx]
    
def reconstruct(old_model, new_model, f):
    macro = calc_macro0(old_model, f)
    return grad13(new_model, macro)

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
        print f[f<0]
        raise NameError("Negative value has been found!")

# Second-order TVD scheme
def transport(domain, bc, solution, delta_t):
    logging.debug('Starting transfer')
    N, xi_y, _xi_y = domain.N, domain.model.xi(np.zeros((domain.N, 3)))[...,1], domain.model.xi0()[...,1]
    gamma, h = _xi_y * delta_t,  delta_y(domain)
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    def calc_F(h, f, idx, F, Idx, mask):
        logging.debug('  - calc_F')
        g = np.einsum('i,a', gamma, 1/h[idx])
        d1, d2 = f[idx+1] - f[idx], f[idx] - f[idx-1]
        h1, h2 = (h[idx+1] + h[idx])/2, (h[idx] + h[idx-1])/2
        # MC limiter
        D = lambda d, h: np.einsum('ai,a->ai', np.abs(d), 1/h)
        H = np.einsum('a,ai->ai', h[idx], np.sign(d1))
        lim = (args.order-1) * H * np.minimum(D(d1+d2, h1+h2), 2*np.minimum(D(d1, h1), D(d2, h2)))
        with np.errstate(divide='ignore', invalid='ignore'):
            F[Idx][mask] = (f[idx] + (1-g)*np.where(d1*d2 > 0, lim, 0)/2)[mask]
    def calc2N(h, f, F, sgn, bc):
        logging.debug(' - calc2N %+d', sgn)
        m0, m1, m2, mN = _mask(sgn), mask(sgn, 1), mask(sgn, 2), mask(sgn, N-2)
        f3[:2][m2], h3[:2] = f[-2:][m2], h[-2:]
        h3[2] = bc.last(h3, f3, f, m0)
        calc_F(h3, f3, np.array([1]), F, slice(-1, None), m1)   # last flux
        bc.last_flux_is_calculated()
        calc_F(h, f, np.arange(1, N-1), F, slice(2, -1), mN)    # interior fluxes
        #check(F[2:][mask(sgn, N-1)]) # after calc2N
    def calc01(h, f, F, sgn, bc):
        logging.debug(' - calc01 %+d', sgn)
        m0, m1, m2 = _mask(sgn), mask(sgn, 1), mask(sgn, 2)
        bc.flux(F, m0, calc_F)                                  # first flux
        f3[1:][m2], h3[1:] = f[:2][m2], h[:2]
        h3[0] = bc.first(h3, f3, f, m0)
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
    f -= np.einsum('i,a,ai->ai', gamma, 1/h, F[1:] - F[:-1])
    check(f) # after transport

def bgk(domain, bc, solution, delta_t):
    logging.debug('Starting BGK')
    macro = calc_macro0(domain.model, solution.f)
    nu = lambda rho: 2/(np.pi**.5*args.kn) * rho
    M = domain.model.Maxw(macro)
    solution.f[:] = np.einsum('ai,a->ai', solution.f - M, np.exp(-delta_t*nu(macro.rho))) + M
    check(solution.f) # after BGK

# Boundary conditions
class Boundary():
    def __init__(self, n):
        pass
    def last_flux_is_calculated(self):
        pass
    def all_fluxes_are_calculated(self):
        pass
    def wait_for_update(self):
        pass
    def last(self, h3, f3, f, mask):            # last ghost cell
        f3[2][mask] = (2*f[-1] - f[-2])[mask]   # extrapolation by default
        return h3[0]
    def flux(self, F, mask, calc_F):            # first flux
        self(F[0], F[0], mask)
    def first(self, h3, f3, f, mask):           # first ghost cell
        self(f3[0], f[0], mask)
        return h3[1]

class Symmetry(Boundary):
    def __call__(self, F, F0, mask):
        F[mask] = F0[::-1][mask]
    def last(self, h3, f3, f, mask):
        self(f3[2], f[-1], mask)
        return h3[1]

class Diffuse(Boundary):
    def __init__(self, n):
        self._model = domains[n].model
        self._xi_y = self._model.xi0()[...,1]
        self._half_M = half_M(self._model)
    def __call__(self, F, F0, mask):
        rho = np.sum((self._xi_y*F0)[::-1][mask])/self._half_M
        F[mask] = self._model.Maxw0(Macro(rho, args.U*e_x/2, fixed.T_B))[mask]

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
        self._f = solution[self._n].f
    def flux(self, F, mask, calc_F):
        model, model_partner = domains[self._n].model, domains[self._n_partner].model
        if model == model_partner:
            # take the prepared flux from the partner
            self._lock_F.wait()
            F[0][mask] = self._F_partner[self._idx_partner][mask]
            self._lock_F.clear()
        else:
            # reconstruct a flux from the partner distribution function (2 cells)
            idx, idx_partner = bc[self._n_partner][self._idx_partner]._idx_partner, self._idx_partner
            f3, h3 = np.empty((3,) + F[0].shape), np.empty(3)
            f3[2], h3[2] = self._f[idx], delta_y(domains[self._n])[idx]
            f3[:2] = reconstruct(model_partner, model, self._f_partner[::(2*idx+1)][-2:])
            h3[:2] = delta_y(domains[self._n_partner])[::(2*idx+1)][-2:]
            calc_F(h3, f3, np.array([1]), F, slice(1), np.array([mask]))
    def first(self, h3, f3, f, mask):
        return self(f3[0], mask)
    def last(self, h3, f3, f, mask):
        return self(f3[2], mask)
    def last_flux_is_calculated(self):
        bc[self._n_partner][self._idx_partner]._lock_F.set()
    def all_fluxes_are_calculated(self):
        bc[self._n_partner][self._idx_partner]._lock_f.set()
    def wait_for_update(self):
        self._lock_f.wait()
        self._lock_f.clear()

def create_bc():
    bc = []
    for n, b in enumerate(boundaries):
        boundary = []
        for p in b:
            kwargs = p[1]
            kwargs['n'] = n
            boundary.append(p[0](**kwargs))
        bc.append(boundary)
    return bc

def new_solution(initial):
    solution = [ Solution(d, initial) for d in domains ]
    for b in bc:
        for boundary in b:
            if boundary.__class__ == Couple:
                boundary.connect(solution)
    return solution

def run_threads(solution, operator, delta_t):
    threads = []
    for d, s, b in zip(domains, solution, bc):
        thread = threading.Thread(target=operator, args=(d, b, s, delta_t))
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
            print '%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel[-1]/U = %g, T[-1] = %g' \
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, m.vel[-1,0]/U, m.temp[-1] )
            plot_profiles(solution)
    total_values(solution)
    y, h, m = calc_macro(solution)
    np.savetxt(sys.stdout, np.transpose((y, m.vel[:,0], m.tau[:,2], m.qflow[:,0])), fmt='%1.4e')
    #splot(domains[-1].model, solution[-1].f[-1])
    #splot(domains[0].model, solution[0].f[0])
    plt.ioff(); plt.show()

def tests():
    def is_almost_equal(a, b):
        assert np.sum(np.abs(a-b)) < 1e-5
    def print_sample(m, m0):
        s = lambda name: m._asdict()[name] - m0._asdict()[name]
        v = lambda name, idx: m._asdict()[name][:,idx] - m0._asdict()[name][idx]
        columns = ( 'rho', 'temp', 'vel_x', 'vel_y', 'p_xy', 'q_x', 'q_y' )
        values = ( s('rho'), s('temp'), v('vel',0), v('vel',1), v('tau',2), v('qflow', 0), v('qflow', 1) )
        print 'Variable: ', ''.join(['%9s ' for i in range(len(columns))]) % columns
        print 'Deviation:', ''.join(['%+.2e ' for i in range(len(values))]) % values
    def check_macro(model, f, expected):
        computed = calc_macro0(model, f)
        print_sample(computed, expected)
        for m, m0 in zip(computed._asdict().values(), expected._asdict().values()):
            is_almost_equal(m, m0)
    o = lambda y: np.ones_like(y)
    delta = args.U
    rho, vel, temp, tau, qflow = 1 + delta, delta*e_x, 1 + delta**2, delta*e_z, delta*e_x/2
    print 'Test #1: Maxwell distribution (rho=%g, vel_x=%g, temp=%g)' % (rho, vel[0], temp)
    macro = Macro(rho, vel, temp)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, model.Maxw0(macro), macro)
    print ''.join(['-' for i in range(50)])
    print 'Test #2: Grad distribution (rho=%g, vel_x=%g, temp=%g, p_xy=%g)' % (rho, vel[0], temp, tau[2])
    macro = Macro(rho, vel, temp, tau)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, grad13(model, macro), macro)
    print ''.join(['-' for i in range(50)])
    print 'Test #3: Grad distribution (rho=%g, vel_x=%g, temp=%g, qflow_x=%g)' % (rho, vel[0], temp, qflow[0])
    macro = Macro(rho, vel, temp, qflow=qflow)
    for name, model in models.iteritems():
        print '-- %s model:' % name
        check_macro(model, grad13(model, macro), macro)

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
    'd3q121': lbm_d3(1.19697977039307435897239, [(0,0,0), (1,0,0), (1,1,1), (1,2,0), (2,2,0), (3,0,0), (2,3,0), (2,2,2), (1,3,0), (3,3,3)],
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
Model = namedtuple('Model', 'info xi0 Maxw0 weights0 xi Maxw')
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
bc = create_bc()

class Solution():
    def __init__(self, domain, initial):
        empty = lambda model, N: np.empty((N,) + model.Maxw0(Macro()).shape)
        self.f = domain.model.Maxw(Macro(*initial(cells(domain))))  # distribution function
        self.F = empty(domain.model, domain.N + 1)                  # fluxes between cells
        self.f3 = empty(domain.model, 3)                            # ghost + 2 cells

print ''.join(['=' for i in range(50)])
print 'DVM: xi_max = %g, grid=(%d)^3, total = %d' % (args.radius, 2*args.M, models['dvm'].xi0().size/3)
print 'LBM: type = %s, total = %d' % (args.lattice, models['lbm'].xi0().size/3)
print 'Kn = %g, U = %g, cells = %d + %d' % (args.kn, args.U, args.N1, args.N2)
print 'Model: (antisym)[ %s | %s ](diffuse), order = %d' % (args.model1, args.model2, args.order)
print 'Width:  |<-- %.3f -->|<-- %.3f -->|' % (fixed.L-args.width, args.width)
print 'Delta_y:     | %.4f | %.4f |' % ((fixed.L-args.width)/args.N1, args.width/args.N2)
print ''.join(['=' for i in range(50)])

if args.tests:
    tests()
else:
    solve_bgk()

