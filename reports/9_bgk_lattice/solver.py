#!/usr/bin/env python
import argparse
import numpy as np
from functools import partial
from collections import namedtuple
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-N', type=int, default=8, metavar='<int>', help='number of cells in the physical space')
parser.add_argument('-M', type=int, default=8, metavar='<int>', help='number of points along each axis in the velocity space')
parser.add_argument('-e', '--end', type=int, default=int(1e2), metavar='<int>', help='maximum total time')
parser.add_argument('-s', '--step', type=float, default=1, metavar='<float>', help='reduce timestep in given times')
parser.add_argument('-r', '--radius', type=float, default=4., metavar='<float>', help='radius of the velocity grid')
parser.add_argument('-m', '--model', default='dvm', metavar='dvm|lbm', help='type of velocity model')
parser.add_argument('-l', '--lattice', default='d3q19', metavar='d3q15|19|27', help='type of velocity lattice')
parser.add_argument('-o', '--order', type=int, default=2, metavar='1|2', help='order of approximation')
parser.add_argument('-k', '--kn', type=float, default=2e-1/np.pi**.5, metavar='<float>', help='the Knudsen number')
parser.add_argument('-U', type=float, default=2e-2, metavar='<float>', help='difference in velocity of plates')
parser.add_argument('-t', '--tests', action='store_true', help='to run some tests')
args = parser.parse_args()

print ''.join(['=' for i in range(20)])
print 'Grid: (%d)x(%d)^3, xi_max = %g' % (args.N, 2*args.M, args.radius)
print 'Kn = %g, U = %g' % (args.kn, args.U)
print 'Model = %s, order = %d' % (args.model, args.order)
print ''.join(['=' for i in range(20)])

### Constants
T_B, L = 1., .5
_zeros, _ones, dom = np.zeros(3), np.ones(3), np.ones(args.N)
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
e_x = np.eye(3)[0]
X = L * (np.arange(args.N) + .5) / args.N               # coordinates of cell centers in physical space

### DVM velocity grid
idx, delta_v = lambda M: np.arange(2*M) - M + .5, args.radius / args.M
_Xi = args.radius / args.M * idx(args.M)
Xi, zeros  = np.einsum('i,a', _Xi, dom), np.einsum('l,a', _zeros, dom)
_I, I = np.ones_like(_Xi), np.ones_like(Xi)

_xi = lambda v: np.einsum('li,lj,lk->ijkl', (_Xi-v[0], _I, _I), (_I, _Xi-v[1], _I), (_I, _I, _Xi-v[2]))
_sqr_xi = lambda v, s=_ones: np.linalg.norm(np.einsum('ijkl,l->ijkl', _xi(v), s), axis=3)**2
_vel = lambda v, i: np.einsum('a,i', v[:,i], _I) 
_Maxw = lambda rho, vel, temp: rho/(np.pi*temp)**1.5 * np.exp(-_sqr_xi(vel)/temp)*delta_v**3

xi = lambda v: np.einsum('lai,laj,lak->aijkl', (Xi-_vel(v,0), I, I), (I, Xi-_vel(v,1), I), (I, I, Xi-_vel(v,2)))
sqr_xi = lambda v: np.einsum('aijkl,aijkl->aijk', xi(v), xi(v))
Maxw = lambda rho, vel, temp: np.einsum('a,aijk->aijk', rho/(np.pi*temp)**1.5*delta_v**3, np.exp(np.einsum('aijk,a->aijk', -sqr_xi(vel), 1./temp)))

half_M = lambda model: -np.sum(model.xi0()[...,1] * model.Maxw0(1., args.U*e_x/2, T_B) * (model.xi0()[...,1] < 0))
_gr1 = lambda vel, temp, tau: 2 * np.einsum('ijkl,ijkm,n,lmn->ijk', xi(vel), xi(vel), tau, hodge)
_gr2 = lambda vel, temp, qflow: 4./5 * (np.einsum('ijkl,ijk,l->ijk', xi(vel), sqr_xi(vel)/temp-2.5, qflow))
Grad13 = lambda rho, vel, temp, tau, qflow: Maxw(vel, temp, rho) * (1 + (_gr1(vel, temp, tau) + _gr2(vel, temp, qflow))/(rho*temp**2))

# 1d packing, symm=1 for antisymmetry, symm=-1 for specular reflection
def dvm_ball(symm):
    _xi_y = _xi(_zeros)[...,1]
    p = np.where((_xi_y > 0)*(_sqr_xi(_zeros) <= args.radius**2))
    n = np.where((_xi_y < 0)*(_sqr_xi(_zeros) <= args.radius**2))
    return np.hstack((n[0], p[0][::symm])), np.hstack((n[1], p[1])), np.hstack((n[2], p[2]))
ball = dvm_ball(1)

### LBM velocity grid
Lattice = namedtuple('Lattice', ['xi', 'w'])
def lbm_d3(nodes, weights):
    lattice, k = [], 0
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
            lattice += [ (node, weights[k]) ]
    for node in nodes:
        _node = np.array(node)
        while True:
            add_symm(0, _node, lattice)
            _node = np.roll(_node, 1)
            if (_node == node).all():
                break
        k += 1
    xi, w = np.transpose(sorted(lattice, key=sort))
    return Lattice(np.vstack(xi), np.hstack(w))

lattices = {
    'd3q15': lbm_d3([(0,0,0), (1,0,0), (1,1,1)], np.array([2, 1, .125])/9),
    'd3q19': lbm_d3([(0,0,0), (1,0,0), (1,1,0)], np.array([3, 18, 36])**-1.),
    'd3q27': lbm_d3([(0,0,0), (1,0,0), (1,1,0), (1,1,1)], np.array([8, 2, .5, .125])/27)
}

_xi_v = lambda v: np.einsum('il,l', lattices[args.lattice].xi, v)
_sqr = lambda v: np.einsum('l,l', v, v)
_weighted = lambda f: np.einsum('i,i->i', lattices[args.lattice].w, f)
xi_v = lambda v: np.einsum('il,al', lattices[args.lattice].xi, v)
sqr = lambda v: np.einsum('al,al,i->ai', v, v, lat)
weighted = lambda rho, f: np.einsum('a,i,ai->ai', rho, lattices[args.lattice].w, f)
lat = np.ones_like(lattices[args.lattice].w)

### Factory of classes
Model = namedtuple('Model', ['xi0', 'Maxw0', 'weights', 'xi', 'Maxw'])
models = {
    'dvm': Model(
        xi0 = lambda v=_zeros: _xi(v)[ball],
        Maxw0 = lambda rho, vel, temp: _Maxw(rho, vel, temp)[ball],
        weights = lambda: np.ones_like(_xy_y) * delta_v**3,
        xi = lambda v=zeros: xi(v)[(slice(None),) + ball],
        Maxw = lambda rho, vel, temp: Maxw(rho, vel, temp)[(slice(None),) + ball]
    ),
    'lbm': Model(
        xi0 = lambda v=_zeros: lattices[args.lattice].xi - np.einsum('i,l', lat, v),
        Maxw0 = lambda rho, vel, temp: rho*_weighted(1 + 3*(_xi_v(vel) + 1.5*_xi_v(vel)**2 - _sqr(vel)/2)),
        weights = lambda: lattices[args.lattice].w,
        xi = lambda v=zeros: np.einsum('il,a', lattices[args.lattice].xi, dom) - np.einsum('i,al', lat, v),
        Maxw = lambda rho, vel, temp: weighted(rho, 1 + 3*(xi_v(vel) + 1.5*xi_v(vel)**2 - sqr(vel)/2))
    )
}

def splot(model, f):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    xi = model.xi0()
    m = xi[:,2] == delta_v/2
    ax.scatter(xi[:,0][m], xi[:,1][m], (f/model.weights())[m])
    plt.show()

def plot_profiles(model, f):
    plt.clf()
    rho, temp, vel, qflow, tau = calc_macro(model, f)
    U = args.U
    x, Vel, Tau = np.loadtxt('k1e-1.txt').T
    plt.plot(X, vel[:,0]/U, 'r', X, -tau[:,2]/U, 'g', x, Vel, 'rD--', x, -2*Tau, 'gD--')
    #plt.plot(X, 1e3*(rho*temp-1), 'b*--')
    plt.plot(x, x, 'r-.', x, 0*x + np.pi**-.5, 'g-.')
    plt.show()
    plt.pause(1e-3)

def calc_macro(model, f):
    rho = np.einsum('ai->a', f)
    vel = np.einsum('ai,il,a->al', f, model.xi0(), 1./rho)
    c, sqr_c = model.xi(vel), np.einsum('ail,ail->ai', model.xi(vel), model.xi(vel))
    csqr_c = np.einsum('ai,ail->ail', sqr_c, c)
    cc_ij = np.einsum('ail,aim->ailm', c, c)
    cc = np.einsum('ailm,lmn', cc_ij, hodge)
    temp = 2./3*np.einsum('ai,ai,a->a', f, sqr_c, 1./rho)
    qflow = np.einsum('ai,ail->al', f, csqr_c)
    tau = 2*np.einsum('ai,ail->al', f, cc)
    return rho, temp, vel, qflow, tau

def total_values(model, f):
    rho, temp, vel, qflow, tau = calc_macro(model, f)
    print 'Total mass =', sum(rho)/args.N

def check(f):
    if np.sum(np.isnan(f)) > 0:
        raise NameError("NaN has been found!")
    if np.sum(f<0) > 0:
        print "Negative: ", f[f<0]
        raise NameError("Negative value has been found!")

def check_maxw(mod, M):
    rho, temp, vel, qflow, tau = calc_macro(mod, M)
    U, fmt1, fmt2 = args.U, '%+.2e ', '%9s '
    columns = ( 'rho-1', 'temp', 'vel_x/U', 'vel_y/U', 'p_xy/U', 'q_x/U' )
    values = ( 0., 1., 1., 0., 0., 0. )
    print 'Values:   ' + ''.join([fmt2 for i in range(len(values))]) % columns
    print 'Expected: ' + ''.join([fmt1 for i in range(len(values))]) % values
    print 'Computed: ' + ''.join([fmt1 for i in range(len(values))]) % \
        ( rho[0]-1, temp[0], 2*vel[-1, 0]/U, 2*vel[-1, 1]/U, tau[-1, 2]/U, qflow[-1, 0]/U )
 
def tests():
    mod = models[args.model]
    check_maxw(mod, mod.Maxw(dom, np.einsum('a,l', dom, args.U*e_x/2), dom))
    check_maxw(mod, np.einsum('a,i', dom, mod.Maxw0(1, args.U*e_x/2, 1)))

# Second-order TVD scheme
def transport(model, f, delta_t, delta_y):
    rho, temp, vel, qflow, tau = calc_macro(model, f)
    _half_M = half_M(model)
    #print "Delta_h = %g, Delta_t = %g" % (delta_y, delta_t)
    N, xi_y, _xi_y = args.N, model.xi()[...,1], model.xi0()[...,1]
    gamma = _xi_y * delta_t / delta_y
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    # Create additional matrices
    shape = np.array(f.shape)
    F = np.zeros(shape + np.array([1, 0]))                      # xi_y*F -- fluxes between cells
    f3 = np.zeros(shape + np.array([3-N, 0]))                   # ghost + two cells
    #F.fill(np.NaN); f3.fill(np.NaN);                            # for debug
    def calc_F(f, idx, F, Idx, mask):
        df1, df2, df3 = [f[idx+1] - f[idx], f[idx] - f[idx-1], f[idx+1] - f[idx-1]]
        # MC limiter
        lim = lambda d1, d2, d3: (args.order-1)*np.minimum(np.abs(d3)/2, 2*np.minimum(np.abs(d1), np.abs(d2)))*np.sign(d1)
        with np.errstate(divide='ignore', invalid='ignore'):
            F[Idx][mask] = (f[idx] + np.einsum('i,ai->ai', (1-gamma)/2, np.where(df1*df2 > 0, lim(df1, df2, df3), 0)))[mask]
    def calc2N(f, F, sgn):
        m0, m1, m2, mN = _mask(sgn), mask(sgn, 1), mask(sgn, 2), mask(sgn, N-2)
        f3[:2][m2] = f[-2:][m2]
        f3[2][m0] = (2*f[-1] - f[-2])[m0]                       # last ghost cell (extrapolation for diffuse)
        if sgn < 0:
            antisym(f3[2], f[-1], m0)                           # last ghost cell (exact for antisymmetry)
        calc_F(f, np.arange(1, N-1), F, slice(2, N), mN)        # interior fluxes
        calc_F(f3, np.arange(1, 2), F, slice(N, N+1), m1)       # last flux
        #check(F[2:][mask(sgn, N-1)])
    def calc01(f, F, sgn, bc):
        m0, m1, m2 = _mask(sgn), mask(sgn, 1), mask(sgn, 2)
        bc(F[0], F[0], m0)                                      # first flux
        f3[1:3][m2] = f[:2][m2]
        bc(f3[0], f[0], m0)
        calc_F(f3, np.arange(1, 2), F, slice(1, 2), m1)         # second flux
        #check(F[:2][m2])
    def antisym(F, F0, mask):
        F[mask] = F0[::-1][mask]
    def diffuse(F, F0, mask):
        F[mask] = model.Maxw0(np.sum((_xi_y*F0)[::-1][mask])/_half_M, args.U*e_x/2, T_B)[mask]
    calc2N(f, F, 1)
    calc2N(f[::-1], F[::-1], -1)
    calc01(f, F, 1, antisym)
    calc01(f[::-1], F[::-1], -1, diffuse)
    #check(f3)
    #check(F)
    f -= np.einsum('i,ai->ai', gamma, F[1:N+1] - F[0:N])
    check(f) # after transport

def bgk(model, f, delta_t):
    rho, temp, vel, qflow, tau = calc_macro(model, f)
    nu = lambda rho: 2/(np.pi**.5*args.kn) * rho
    M = model.Maxw(rho, vel, temp)
    f[:] = np.einsum('ai,a->ai', f - M, np.exp(-delta_t*nu(rho))) + M
    check(f) # after BGK

def solve_bgk(model):
    delta_y, U, N = L/args.N, args.U, args.N
    U0 = .9*U if args.kn < 1 else 0
    Rho0, Vel0, Temp0 = dom, np.einsum('a,l', X, U0*e_x), dom
    f = model.Maxw(Rho0, Vel0, Temp0)
    rho0, temp0, vel0, qflow0, tau0 = calc_macro(model, f)
    total_values(model, f)
    plt.ion()
    plot_profiles(model, f)
    delta_t = delta_y/args.radius/args.step                     # from CFL
    for i in xrange(args.end):
        # Symmetry second-order splitting
        transport(model, f, delta_t, delta_y)
        bgk(model, f, 2*delta_t)
        transport(model, f, delta_t, delta_y)
        if not i % 10:
            rho, temp, vel, qflow, tau = calc_macro(model, f)
            rho_disp, pxy_mean = sum((rho-rho0)**2)/N, np.sum(tau[:,2])/N
            pxy_disp = sum((tau[:,2]-pxy_mean)**2)/N
            print '%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel[-1]/U = %g, T[-1] = %g' \
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, vel[-1,0]/U, temp[-1] )
            plot_profiles(model, f)
    total_values(model, f)
    rho, temp, vel, qflow, tau = calc_macro(model, f)
    print rho, temp, vel/U, qflow/U, tau/U
    #splot(f[-1])
    #splot(f[0])
    plt.ioff(); plt.show()

if args.tests:
    tests()
else:
    solve_bgk(models[args.model])

