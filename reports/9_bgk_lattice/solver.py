#!/usr/bin/env python
import argparse, threading, logging
import numpy as np
from functools import partial
from collections import namedtuple
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-N1', type=int, default=4, metavar='<int>', help='number of cells in the first domain')
parser.add_argument('-N2', type=int, default=4, metavar='<int>', help='number of cells in the second domain')
parser.add_argument('-M', type=int, default=8, metavar='<int>', help='number of points along each axis in the velocity space')
parser.add_argument('-e', '--end', type=int, default=int(1e2), metavar='<int>', help='maximum total time')
parser.add_argument('-w', '--width', type=float, default=1.2, metavar='<float>', help='width of the second domain in terms of mean free path')
parser.add_argument('-s', '--step', type=float, default=1, metavar='<float>', help='reduce timestep in given times')
parser.add_argument('-r', '--radius', type=float, default=4., metavar='<float>', help='radius of the velocity grid')
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

if args.verbose:
    logging.basicConfig(level=logging.DEBUG, format='(%(threadName)-10s) %(message)s')

### Constants
T_B, L = 1., .5
width = min(args.width*args.kn, L/2)
zeros, ones = np.zeros(3), np.ones(3)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
e_x = np.eye(3)[0]

### Auxiliary functions
cells = lambda domain: domain.y0 + domain.L*(np.arange(domain.N) + .5)/domain.N
delta_y = lambda domain: domain.L/domain.N
half_M = lambda model: -np.sum(model.xi0()[...,1] * model.Maxw0(1., args.U*e_x/2, T_B) * (model.xi0()[...,1] < 0))

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
    _Maxw = lambda rho, vel, temp: rho/(np.pi*temp)**1.5 * np.exp(-_sqr_xi(vel)/temp)*delta_v**3

    xi = lambda v: np.einsum('lai,laj,lak->aijkl', G(v, 0), G(v, 1), G(v, 2))
    sqr_xi = lambda v: np.einsum('aijkl,aijkl->aijk', xi(v), xi(v))
    Maxw = lambda rho, vel, temp: np.einsum('a,aijk->aijk', rho/(np.pi*temp)**1.5*delta_v**3, np.exp(np.einsum('aijk,a->aijk', -sqr_xi(vel), 1./temp)))

    # 1d packing, symm=1 for antisymmetry, symm=-1 for specular reflection
    def dvm_ball(symm):
        _xi_y = _xi(zeros)[...,1]
        p = np.where((_xi_y > 0)*(_sqr_xi(zeros) <= args.radius**2))
        n = np.where((_xi_y < 0)*(_sqr_xi(zeros) <= args.radius**2))
        return np.hstack((n[0], p[0][::symm])), np.hstack((n[1], p[1])), np.hstack((n[2], p[2]))
    ball = dvm_ball(1)

    return Model(
        xi0 = lambda v=zeros: _xi(v)[ball],
        Maxw0 = lambda rho=1, vel=zeros, temp=1: _Maxw(rho, vel, temp)[ball],
        weights0 = lambda: np.ones_like(ball[0]) * delta_v**3,
        xi = lambda v: xi(v)[(slice(None),) + ball],
        Maxw = lambda rho, vel, temp: Maxw(rho, vel, temp)[(slice(None),) + ball]
    )

def grad13():
    _gr1 = lambda vel, temp, tau: 2 * np.einsum('ijkl,ijkm,n,lmn->ijk', xi(vel), xi(vel), tau, hodge)
    _gr2 = lambda vel, temp, qflow: 4./5 * (np.einsum('ijkl,ijk,l->ijk', xi(vel), sqr_xi(vel)/temp-2.5, qflow))
    Grad13 = lambda rho, vel, temp, tau, qflow: Maxw(vel, temp, rho) * (1 + (_gr1(vel, temp, tau) + _gr2(vel, temp, qflow))/(rho*temp**2))

def lbm_d3(nodes, weights):
    lattice = []
    def sort(item):
        return 2*np.sum(item[0]) + np.sum(item[0]/2**np.arange(3.))
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
    for node in nodes:
        _node = np.array(node)
        while True:
            add_symm(0, _node, lattice)
            _node = np.roll(_node, 1)
            if (_node == node).all():
                break
        weights = np.roll(weights, -1)
    xi, w = np.transpose(sorted(lattice, key=sort))
    return Lattice(np.vstack(xi), np.hstack(w))

def lbm_grid():
    _xi_v = lambda v: np.einsum('il,l', lattices[args.lattice].xi, v)
    _sqr = lambda v: np.einsum('l,l', v, v)
    _weighted = lambda f: np.einsum('i,i->i', lattices[args.lattice].w, f)
    xi_v = lambda v: np.einsum('il,al', lattices[args.lattice].xi, v)
    sqr = lambda v: np.einsum('al,al,i->ai', v, v, lat)
    weighted = lambda rho, f: np.einsum('a,i,ai->ai', rho, lattices[args.lattice].w, f)
    lat = np.ones_like(lattices[args.lattice].w)
    Maxw0 = lambda rho=1, vel=zeros, temp=1: rho*_weighted(1 + 3*(_xi_v(vel) + 1.5*_xi_v(vel)**2 - _sqr(vel)/2))
    return Model(
        xi0 = lambda v=zeros: lattices[args.lattice].xi - np.einsum('i,l', lat, v),
        Maxw0 = lambda rho=1, vel=zeros, temp=1: rho*_weighted(1 + 3*(_xi_v(vel) + 1.5*_xi_v(vel)**2 - _sqr(vel)/2)),
        weights0 = lambda: lattices[args.lattice].w,
        xi = lambda v: np.einsum('il,a', lattices[args.lattice].xi, np.ones(v.shape[0])) - np.einsum('i,al', lat, v),
        Maxw = lambda rho, vel, temp: weighted(rho, 1 + 3*(xi_v(vel) + 1.5*xi_v(vel)**2 - sqr(vel)/2))
    )

def splot(model, f):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    xi = model.xi0()
    m = ( xi[:,2] <= (model.weights0()[0])**(1./3)/2 ) * ( xi[:,2] >= 0 ) 
    #ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.weights0())[m])
    ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.Maxw0())[m])
    plt.show()

def plot_profiles(solution):
    plt.clf()
    y, h, rho, vel, temp, tau, qflow = calc_macro(solution)
    U = args.U
    Y, Vel, Tau = np.loadtxt('k1e-1.txt').T
    plt.plot(Y, Vel, 'rD--', Y, -2*Tau, 'gD--')             # k = 0.1
    Y, Vel, Tau = np.loadtxt('k1e-0.txt').T
    plt.plot(Y, Vel, 'rs--', Y, -2*Tau, 'gs--')             # k = 1
    plt.plot(Y, Y, 'r-.')                                   # k = 0
    plt.plot(Y, 0*Y + np.pi**-.5, 'g-.')                    # k = \infty
    plt.plot(y, vel[:,0]/U, 'r*-', y, -tau[:,2]/U, 'g*-')   # present
    plt.show()
    plt.pause(1e-3)

def calc_macro(solution):
    y, h, rho, vel, temp, tau, qflow = [], [], [], [], [], [], []
    #TODO: optimize using generators
    for d, s in zip(domains, solution):
        y += [ cells(d) ]
        h += [ delta_y(d)*np.ones(d.N) ]
        rho += [ np.einsum('ai->a', s.f) ]
        vel += [ np.einsum('ai,il,a->al', s.f, d.model.xi0(), 1./rho[-1]) ]
        c = d.model.xi(vel[-1])
        sqr_c = np.einsum('ail,ail->ai', c, c)
        csqr_c = np.einsum('ai,ail->ail', sqr_c, c)
        cc_ij = np.einsum('ail,aim->ailm', c, c)
        cc = np.einsum('ailm,lmn', cc_ij, hodge)
        temp += [ 2./3*np.einsum('ai,ai,a->a', s.f, sqr_c, 1./rho[-1]) ]
        qflow += [ np.einsum('ai,ail->al', s.f, csqr_c) ]
        tau += [ 2*np.einsum('ai,ail->al', s.f, cc) ]
    s = lambda m: np.hstack(m)
    v = lambda m: np.vstack(m)
    return s(y), s(h), s(rho), v(vel), s(temp), v(tau), v(qflow)

def total_values(domains):
    y, h, rho, vel, temp, tau, qflow = calc_macro(domains)
    print 'Total mass =', 2*sum(h*rho)

def check(f):
    if np.sum(np.isnan(f)) > 0:
        raise NameError("NaN has been found!")
    if np.sum(f<0) > 0:
        print "Negative: ", f[f<0]
        raise NameError("Negative value has been found!")

# Second-order TVD scheme
def transport(domain, bc, solution, delta_t):
    logging.debug('Starting transfer')
    N, xi_y, _xi_y = domain.N, domain.model.xi(np.zeros((domain.N, 3)))[...,1], domain.model.xi0()[...,1]
    gamma = _xi_y * delta_t / delta_y(domain)
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    #F.fill(np.NaN); f3.fill(np.NaN);                            # for debug
    def calc_F(f, idx, F, Idx, mask):
        logging.debug('  - calc_F')
        df1, df2, df3 = [f[idx+1] - f[idx], f[idx] - f[idx-1], f[idx+1] - f[idx-1]]
        # MC limiter
        lim = lambda d1, d2, d3: (args.order-1)*np.minimum(np.abs(d3)/2, 2*np.minimum(np.abs(d1), np.abs(d2)))*np.sign(d1)
        with np.errstate(divide='ignore', invalid='ignore'):
            F[Idx][mask] = (f[idx] + np.einsum('i,ai->ai', (1-gamma)/2, np.where(df1*df2 > 0, lim(df1, df2, df3), 0)))[mask]
    def calc2N(f, F, sgn, bc):
        logging.debug(' - calc2N')
        m0, m1, m2, mN = _mask(sgn), mask(sgn, 1), mask(sgn, 2), mask(sgn, N-2)
        f3[:2][m2] = f[-2:][m2]
        f3[2][m0] = (2*f[-1] - f[-2])[m0]                       # last ghost cell (extrapolation for diffuse)
        if bc.__class__ == Symmetry:
            bc.cell(f3[2], f[-1], m0)                           # last ghost cell (exact for antisymmetry)
        calc_F(f, np.arange(1, N-1), F, slice(2, N), mN)        # interior fluxes
        calc_F(f3, np.arange(1, 2), F, slice(N, N+1), m1)       # last flux
        bc.update()
        #check(F[2:][mask(sgn, N-1)])
    def calc01(f, F, sgn, bc):
        logging.debug(' - calc01')
        m0, m1, m2 = _mask(sgn), mask(sgn, 1), mask(sgn, 2)
        bc.flux(F[0], F[0], m0)                                 # first flux
        f3[1:3][m2] = f[:2][m2]
        bc.cell(f3[0], f[0], m0)
        calc_F(f3, np.arange(1, 2), F, slice(1, 2), m1)         # second flux
        #check(F[:2][m2])
    f, F, f3 = solution.f, solution.F, solution.f3
    calc2N(f, F, 1, bc[1])
    calc2N(f[::-1], F[::-1], -1, bc[0])
    calc01(f, F, 1, bc[0])
    calc01(f[::-1], F[::-1], -1, bc[1])
    #check(f3)
    #check(F)
    f -= np.einsum('i,ai->ai', gamma, F[1:N+1] - F[0:N])
    check(f) # after transport

def bgk(domain, bc, solution, delta_t):
    logging.debug('Starting BGK')
    y, h, rho, vel, temp, tau, qflow = calc_macro([solution])
    nu = lambda rho: 2/(np.pi**.5*args.kn) * rho
    M = domain.model.Maxw(rho, vel, temp)
    solution.f[:] = np.einsum('ai,a->ai', solution.f - M, np.exp(-delta_t*nu(rho))) + M
    check(solution.f) # after BGK

# Boundary conditions
class Boundary():
    def __init__(self, n):
        pass
    def update(self):
        pass
    def flux(self, F, F0, mask):
        self(F, F0, mask)
    def cell(self, f, f0, mask):
        self(f, f0, mask)

class Symmetry(Boundary):
    def __call__(self, F, F0, mask):
        F[mask] = F0[::-1][mask]

class Diffuse(Boundary):
    def __init__(self, n):
        self._model = domains[n].model
        self._xi_y = self._model.xi0()[...,1]
        self._half_M = half_M(self._model)
    def __call__(self, F, F0, mask):
        rho = np.sum((self._xi_y*F0)[::-1][mask])/self._half_M
        F[mask] = self._model.Maxw0(rho, args.U*e_x/2, T_B)[mask]

class Couple(Boundary):
    def __init__(self, n, n_partner, idx):
        self._n_partner = n_partner
        self._idx = idx
        self._coeff = delta_y(domains[n]) / delta_y(domains[n_partner])
        self._lock = threading.Event()
    def connect(self, solution):
        self._lock_partner = bc[self._n_partner][self._idx]._lock
        self._f = solution[self._n_partner].f
        self._F = solution[self._n_partner].F
    def flux(self, F, F0, mask):
        self._lock.wait()
        F[mask] = self._F[self._idx][mask]
        self._lock.clear()
    def cell(self, f, f0, mask):
        f[mask] = (self._f[self._idx]*self._coeff + f0*(1-self._coeff))[mask]
    def update(self):
        self._lock_partner.set()

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
        thread = threading.Thread(target=operator, args=(d, b, s, delta_t), name=operator.__name__)
        thread.start()
        threads.append(thread)
    for t in threads:
        t.join()            # wait until all threads terminate

def solve_bgk():
    delta_t = min(map(delta_y, domains))/args.radius/args.step
    U, U0 = args.U, .9*args.U if args.kn < 1 else 0
    solution = new_solution(lambda y: (1+0*y, U0*np.einsum('a,l', y, e_x), 1+0*y))
    y, h, rho0, vel0, temp0, tau0, qflow0 = calc_macro(solution)
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
            y, h, rho, vel, temp, tau, qflow = calc_macro(solution)
            rho_disp, pxy_mean = sum(h*(rho-rho0)**2), np.sum(h*tau[:,2])
            pxy_disp = sum(h*(tau[:,2]-pxy_mean)**2)
            print '%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel[-1]/U = %g, T[-1] = %g' \
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, vel[-1,0]/U, temp[-1] )
            plot_profiles(solution)
    total_values(solution)
    y, h, rho, vel, temp, tau, qflow = calc_macro(solution)
    #print rho, temp, vel/U, qflow/U, tau/U
    #splot(domains[-1].model, solution[-1].f[-1])
    #splot(domains[0].model, solution[0].f[0])
    plt.ioff(); plt.show()

def check_maxw(solution):
    y, h, rho, vel, temp, tau, qflow = calc_macro(solution)
    U, fmt1, fmt2 = args.U, '%+.2e ', '%9s '
    columns = ( 'rho-1', 'temp', 'vel_x/U', 'vel_y/U', 'p_xy/U', 'q_x/U' )
    values = ( 0., 1., 1., 0., 0., 0. )
    print 'Values:   ' + ''.join([fmt2 for i in range(len(values))]) % columns
    print 'Expected: ' + ''.join([fmt1 for i in range(len(values))]) % values
    print 'Computed: ' + ''.join([fmt1 for i in range(len(values))]) % \
        ( rho[0]-1, temp[0], 2*vel[-1, 0]/U, 2*vel[-1, 1]/U, tau[-1, 2]/U, qflow[-1, 0]/U )

def tests():
    o = lambda y: np.ones_like(y)
    solution = new_solution(lambda y: (o(y), args.U/2*np.einsum('a,l', o(y), e_x), o(y)))
    check_maxw(solution)

Lattice = namedtuple('Lattice', ['xi', 'w'])
lattices = {
    'd3q15': lbm_d3([(0,0,0), (1,0,0), (1,1,1)], np.array([2, 1, .125])/9),
    'd3q19': lbm_d3([(0,0,0), (1,0,0), (1,1,0)], np.array([3, 18, 36])**-1.),
    'd3q27': lbm_d3([(0,0,0), (1,0,0), (1,1,0), (1,1,1)], np.array([8, 2, .5, .125])/27)
}

Model = namedtuple('Model', ['xi0', 'Maxw0', 'weights0', 'xi', 'Maxw'])
models = {
    'dvm': dvm_grid(),
    'lbm': lbm_grid()
}

Domain = namedtuple('Domain', ['y0', 'L', 'model', 'N' ])
domains = (
    Domain(0, L-width, models[args.model1], args.N1),
    Domain(L-width, width, models[args.model2], args.N2)
)
boundaries = (
    [ ( Symmetry, {} ),  ( Couple, { 'n_partner': 1, 'idx': 0 }) ],
    [ ( Couple, { 'n_partner': 0, 'idx': -1 }), ( Diffuse, {}) ]
)
bc = create_bc()

class Solution():
    def __init__(self, domain, initial):
        empty = lambda model, N: np.empty((N,) + model.Maxw0().shape)
        self.f = domain.model.Maxw(*initial(cells(domain)))
        self.F = empty(domain.model, domain.N + 1)                  # fluxes between cells
        self.f3 = empty(domain.model, 3)                            # ghost + 2 cells

print ''.join(['=' for i in range(50)])
print 'DVM: xi_max = %g, grid=(%d)^3, total = %d' % (args.radius, 2*args.M, models['dvm'].xi0().size)
print 'LBM: type = %s' % (args.lattice)
print 'Kn = %g, U = %g, cells = %d + %d' % (args.kn, args.U, args.N1, args.N2)
print 'Model: (antisym)[ %s | %s ](diffuse), order = %d' % (args.model1, args.model2, args.order)
print 'Width:  |<-- %.3f -->|<-- %.3f -->|' % (L-width, width)
print 'Delta_y:     | %.4f | %.4f |' % ((L-width)/args.N1, width/args.N2)
print ''.join(['=' for i in range(50)])

if args.tests:
    tests()
else:
    solve_bgk()

