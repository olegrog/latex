#!/usr/bin/env python
r'Solver for the plane Couette-flow problem'
r'Example: ./solver.py --solver=bgk -N=10 -M=10'

import argparse
import numpy as np
import pylab as py

parser = argparse.ArgumentParser(description='2D contour plots of scalar and vector from VTK data')
parser.add_argument('-N', type=int, default=8, metavar='value', help='number of cells in the physical space')
parser.add_argument('-M', type=int, default=8, metavar='value', help='number of points along each axis in the velocity space')
parser.add_argument('-t', '--time', type=int, default=int(1e2), metavar='value', help='Maximum total time')
parser.add_argument('--step', type=float, default=1, metavar='value', help='Reduce timestep in given times')
parser.add_argument('-r', '--radius', type=float, default=4., metavar='value', help='radius of the velocity grid')
parser.add_argument('-s', '--solver', default='bgk', metavar='value', help='type of solver')
parser.add_argument('--scheme', default='implicit', metavar='value', help='explicit or implicit')
parser.add_argument('--order', type=int, default=2, metavar='1|2', help='order of approximation')
parser.add_argument('-k', '--kn', type=float, default=5e-1*np.pi**.5/2, metavar='value', help='Knudsen number')
parser.add_argument('-U', type=float, default=4e-2, metavar='value', help='difference in velocity of plates')
args = parser.parse_args()

print ''.join(['=' for i in range(20)])
print 'Grid: (%d)x(%d)^3, xi_max = %g' % (args.N, 2*args.M, args.radius)
print 'Kn = %g, U = %g' % (args.kn, args.U)
print 'Solver = %s, scheme = %s, order = %d' % (args.solver, args.scheme, args.order)
print ''.join(['=' for i in range(20)])

### Constants
zeros, ones, dom, L = np.zeros(3), np.ones(3), np.ones(args.N), .5
hodge = np.zeros((3,3,3))
hodge[0, 1, 2] = hodge[1, 2, 0] = hodge[2, 0, 1] = 1
e_x, e_y, e_z = np.array([1.,0,0]), np.array([0,1.,0]), np.array([0,0,1.])
X = L * (np.arange(args.N) + .5) / args.N

### Velocity grid
idx = lambda M: np.arange(2*M) - M + .5
_H, _Xi = args.radius / args.M + idx(args.M)*0, args.radius / args.M * idx(args.M)
H, Xi = np.einsum('i,a', _H, dom), np.einsum('i,a', _Xi, dom)
_I, I, Zeros = np.ones_like(_Xi), np.ones_like(Xi), np.einsum('l,a', zeros, dom)

dxi = np.einsum('i,j,k', _H, _H, _H)
_xi = lambda v=zeros: np.einsum('li,lj,lk->ijkl', (_Xi-v[0], _I, _I), (_I, _Xi-v[1], _I), (_I, _I, _Xi-v[2]))
_sqr_xi = lambda vel, s=ones: np.linalg.norm(np.einsum('ijkl,l->ijkl', _xi(vel), s), axis=3)**2
_vel = lambda v, i: np.einsum('a,i', v[:,i], _I) 
_Maxwell = lambda rho, vel, temp: rho/(np.pi*temp)**1.5 * np.exp(-_sqr_xi(vel)/temp)*dxi

xi = lambda v=Zeros: np.einsum('lai,laj,lak->aijkl', (Xi-_vel(v,0), I, I), (I, Xi-_vel(v,1), I), (I, I, Xi-_vel(v,2)))
sqr_xi = lambda vel, s=ones: np.linalg.norm(np.einsum('aijkl,l->aijkl', xi(vel), s), axis=4)**2
Maxwell = lambda rho, vel, temp: np.einsum('a,aijk,ijk->aijk', rho/(np.pi*temp)**1.5, np.exp(np.einsum('aijk,a->aijk', -sqr_xi(vel), 1./temp)), dxi)

_gr1 = lambda vel, temp, tau: 2 * np.einsum('ijkl,ijkm,n,lmn->ijk', xi(vel), xi(vel), tau, hodge)
_gr2 = lambda vel, temp, qflow: 4./5 * (np.einsum('ijkl,ijk,l->ijk', xi(vel), sqr_xi(vel)/temp-2.5, qflow))
Grad13 = lambda rho, vel, temp, tau, qflow: Maxwell(vel, temp, rho) * (1 + (_gr1(vel, temp, tau) + _gr2(vel, temp, qflow))/(rho*temp**2))

# Diffuse boundary constants
T_B = 1.
half_M = -np.sum(_xi()[...,1] * _Maxwell(1., args.U*e_x/2, T_B) * (_xi()[...,1] < 0))

def splot(f):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = py.figure()
    ax = p3.Axes3D(fig)
    xv, yv = np.meshgrid(_Xi, _Xi, sparse=False, indexing='ij')
    f /= dxi
    ball = _sqr_xi(zeros, ones/args.radius) <= 1
    f[np.invert(ball)] = np.NaN
    ax.plot_wireframe(xv, yv, f[:,:,args.M], linewidth=0.25)
    py.show()

def plot_profiles(f):
    py.clf()
    rho, temp, vel, qflow, tau = calc_macro(f)
    U = args.U
    x, Vel, Tau = np.loadtxt('k1e-1.txt').T
    py.plot(X, vel[:,0]/U, 'r', X, -tau[:,2]/U, 'g', x, Vel, 'rD--', x, -Tau, 'gD--')
    #py.plot(X, 1e3*(rho*temp-1), 'b*--')
    py.plot(x, x, 'r-.', x, 0*x + np.pi**-.5, 'g-.')
    py.show()
    py.pause(1e-3)

def calc_macro(f):
    #f = np.einsum('aijk,ijk->ijk', f, ball)
    rho = np.einsum('aijk->a', f)
    vel = np.einsum('aijk,ijkl,a->al', f, _xi(), 1./rho)
    c, sqr_c = xi(vel), sqr_xi(vel)
    csqr_c = np.einsum('aijk,aijkl->aijkl', sqr_c, c)
    cc_ij = np.einsum('aijkl,aijkm->aijklm', c, c)
    cc = np.einsum('aijklm,lmn', cc_ij, hodge)
    temp = 2./3*np.einsum('aijk,aijk,a->a', f, sqr_c, 1./rho)
    qflow = np.einsum('aijk,aijkl->al', f, csqr_c)
    tau = 2*np.einsum('aijk,aijkl->al', f, cc)
    return rho, temp, vel, qflow, tau

def total_values(f):
    rho, temp, vel, qflow, tau = calc_macro(f)
    print 'Total mass =', sum(rho)/args.N

def test():
    sca_ones, vec_zeros = np.array([1.,]), np.array([zeros,])
    rho, temp, vel, qflow, tau = calc_macro(Maxwell(sca_ones, vec_zeros, sca_ones))
    print rho[0], temp[0], vel[0], qflow[0], tau[0]
    # print calc_macro(Grad13([0,0,1e-1],zeros,[1e-1,0,0],1,1))

def check(f):
    if np.sum(np.isnan(f)) > 0:
        raise NameError("NaN has been found!")
    '''
    eps = 1e-4
    f[(f<0)*(f>-eps)] = 0
    if np.sum(f<0) > 0:
        print "Negative: ", f[f<=0]
        raise NameError("Negative value has been found!")
    '''

# Second-order TVD scheme
def iterate_bgk(f, delta_y):
    nu = lambda rho: 2/(np.pi**.5*args.kn) * rho
    check(f)
    rho, temp, vel, qflow, tau = calc_macro(f)
    delta_t = delta_y/args.radius/args.order/2/args.step        # from CFL
    #print "Delta_h = %g, Delta_t = %g" % (delta_y, delta_t)
    N, xi_y, _xi_y = args.N, xi()[...,1], _xi()[...,1]
    gamma = _xi_y * delta_t / delta_y
    mask, _mask = lambda sgn, N: sgn*xi_y[0:N] > 0, lambda sgn: sgn*_xi_y > 0
    # Create additional matrices
    shape = np.array(f.shape)
    F = np.empty(shape + np.array([1, 0, 0, 0]))                # xi_y*F -- fluxes between cells
    f3 = np.empty(shape + np.array([3-N, 0, 0, 0]))             # ghost + two cells
    F.fill(np.NaN); f3.fill(np.NaN);                            # for debug
    def calc_F(f, idx, F, Idx, mask):
        df1, df2, df3 = [f[idx+1] - f[idx], f[idx] - f[idx-1], f[idx+1] - f[idx-1]]
        check(df1[mask])
        check(df2[mask])
        check(df3[mask])
        # MC limiter
        lim = lambda d1, d2, d3: (args.order-1)*np.minimum(np.abs(d3)/2, 2*np.minimum(np.abs(d1), np.abs(d2)))*np.sign(d1)
        F[Idx][mask] = (f[idx] + np.einsum('ijk,aijk->aijk', (1-gamma)/2, np.where(df1*df2 > 0, lim(df1, df2, df3), 0)))[mask]
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
        #f3[0][m0] = (2*f3[0] - f[0])[m0]
        calc_F(f3, np.arange(1, 2), F, slice(1, 2), m1)         # second flux
        #check(F[:2][m2])
    def antisym(F, F0, mask):
        F[mask] = F0[::-1,::-1][mask]
    def diffuse(F, F0, mask):
        F[mask] = _Maxwell(np.sum((_xi_y*F0)[:,::-1][mask])/half_M, args.U*e_x/2, T_B)[mask]
    calc2N(f, F, 1)
    calc2N(f[::-1], F[::-1], -1)
    calc01(f, F, 1, antisym)
    calc01(f[::-1], F[::-1], -1, diffuse)
    check(f3)
    check(F)
    nu_f = lambda rho, vel, temp, f: delta_t*np.einsum('a, aijk', nu(rho), Maxwell(rho, vel, temp) - f)
    fluxes = np.einsum('ijk,aijk->aijk', gamma, F[1:N+1] - F[0:N])
    if args.scheme == 'implicit':
        rho, temp, vel, qflow, tau = calc_macro(f-fluxes)
        #print rho, temp, vel
        vel[:,1] = vel[:,2] = 0
        f += np.einsum('aijk,a->aijk',  nu_f(rho, vel, temp, f) - fluxes, 1./(1 + delta_t*nu(rho)))
    else:
        f += nu_f(rho, vel, temp, f) - fluxes

def solve_bgk():
    delta_y, U, N = L/args.N, args.U, args.N
    U0 = .9*U if args.kn < 1 else 0
    Rho0, Vel0, Temp0 = dom, np.einsum('a,l',X, U0*e_x), dom
    f = Maxwell(Rho0, Vel0, Temp0)
    rho0, temp0, vel0, qflow0, tau0 = calc_macro(f)
    total_values(f)
    py.ion()
    plot_profiles(f)
    for i in xrange(args.time):
        iterate_bgk(f, delta_y)
        if not i % 10:
            rho, temp, vel, qflow, tau = calc_macro(f)
            rho_disp, pxy_mean = sum((rho-rho0)**2)/N, np.sum(tau[:,2])/N
            pxy_disp = sum((tau[:,2]-pxy_mean)**2)/N
            print '%d: err(rho) = %g, p_xy/U = %g, err(p_xy)/U = %g, vel[-1]/U = %g, T[-1] = %g' \
                % ( i, rho_disp, pxy_mean/U, pxy_disp/U, vel[-1,0]/U, temp[-1] )
            plot_profiles(f)
    total_values(f)
    rho, temp, vel, qflow, tau = calc_macro(f)
    print rho, temp, vel/U, qflow/U, tau/U
    #splot(f[-1])
    splot(f[0])
    py.ioff(); py.show()

{
    'bgk': solve_bgk,
    'lbe': solve_bgk,
    'hybrid': solve_bgk,
}[args.solver]()
