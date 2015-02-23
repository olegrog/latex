#!/usr/bin/python

import pylab as py
import numpy as np
import numpy.random
import sys
from scipy.special import erf
from scipy.interpolate import interp1d
from scipy.integrate import quad
from functools import partial
from scipy.optimize import curve_fit

a1, a2, a3 = 1./(np.sqrt(2)*np.pi), np.sqrt(np.pi)/2, 8**(-.5)
N_R, N_V, cut = 27*2+1, 5e4, 5.4
if (len(sys.argv) > 3):
    N_V = int(float(sys.argv[3]))
V_I = 4./3 * np.pi * cut**3
X = np.linspace(0, cut, N_R)

gamma_1 = 1.270042427
sqr = lambda x: np.sum(np.square(x), 0) 
mag = lambda x: np.sqrt(sqr(x))
zeta1d = lambda x: x*np.array([1.,0.,0.])
zeta2d = lambda x: x*np.array([1.,1.,0.])/2**.5
zeta3d = lambda x: x*np.array([1.,1.,1.])/3**.5 # * 0.99999999999

def L1(phi, zeta, xi):
    return a1/mag(zeta-xi) * np.exp(sqr(np.cross(zeta,xi,axis=0))/sqr(zeta-xi) - sqr(xi)) * phi(xi)

def L2(phi, zeta, xi):
    return .5*a1*mag(zeta-xi) * np.exp(-sqr(xi)) * phi(xi)

def nu(Zeta):
    return 2**(-1.5) * (2 if Zeta==0 else (np.exp(-np.square(Zeta)) + (2*Zeta+1./Zeta)*erf(Zeta)*a2))

def get_grid1(N):
    return np.random.rand(N, 3)

def show():
    sys.stdout.flush()
    sys.stderr.flush()
    py.xlim(0, cut)
    py.show()

korobov3 = [
    [ 10007,    [ 1, 544,     5733    ]],
    [ 50021,    [ 1, 12962,   42926   ]],
    [ 100003,   [ 1, 47283,   15021   ]],
    [ 500009,   [ 1, 33606,   342914  ]],
    [ 1000003,  [ 1, 342972,  439897  ]],
    [ 2000003,  [ 1, 235672,  1208274 ]],
    [ 10000019, [ 1, 4341869, 594760  ]]
]

korobov5 = [
    [ 100003,   [ 1, 11729,   65316,   68384,   51876   ]],
    [ 200003,   [ 1, 62638,   60193,   112581,  142904  ]],
    [ 500009,   [ 1, 191775,  488648,  283447,  69999   ]],
    [ 1000003,  [ 1, 335440,  656043,  403734,  126676  ]],
    [ 2000003,  [ 1, 701679,  680513,  965077,  1248525 ]],
    [ 10000019, [ 1, 3669402, 5455092, 7462912, 2188321 ]],
]

def get_grid2(N, dim=3):
    korobov = {
        3: korobov3,
        5: korobov5
    }[dim]
    N, coeffs = [line for line in korobov if line[0] > N][0]
    return np.modf(np.outer(1.+np.arange(N), coeffs) / N + np.random.rand(dim))[0]

def get_xi(N):
    xi = 2*get_grid2(N) - 1
    return (cut*xi[np.sum(np.square(xi), axis=1) < 1]).T

def get_xi_alpha():
    grid = get_grid2(N_V, 5)
    xi = 2*grid[:,0:3] - 1
    mask = np.sum(np.square(xi), axis=1) < 1
    theta, phi = np.pi*grid[:,3][mask], 2*np.pi*grid[:,4][mask]
    s0, c0 = np.sin(theta), np.cos(theta)
    s1, c1 = np.sin(phi), np.cos(phi)
    return (cut*xi[mask]).T, np.array([s0*c1,s0*s1,c0]), s0

def splot(X, Y, K):
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = py.figure()
    ax = p3.Axes3D(fig)
    f = np.abs(K) > 2e-2
    ax.scatter(X[f],Y[f],K[f],s=5,c=K[f])
#    xv, yv = np.meshgrid(X, X, sparse=False, indexing='ij')
#    ax.plot_wireframe(xv, yv, K, rstride=N/10, cstride=N/10)
    show()

def eval_L3(phi, zeta):
    N = N_V
#    if mag(zeta) < .8:
#        N = 10*N_V
    xi = get_xi(N)
    N_I = xi.shape[1]
    #print N_I
    tile_zeta = np.tile(zeta, (N_I, 1)).T
    return (L1(phi, tile_zeta, xi) - L2(phi, tile_zeta, xi)).sum()*V_I/N_I - nu(mag(zeta))*phi(zeta)
    
def eval_L_plus_nu(phi, zeta):
    N = N_V
#    if mag(zeta) < .8:
#        N = 10*N_V
    xi = get_xi(N)
    N_I = xi.shape[1]
    #print N_I
    tile_zeta = np.tile(zeta, (N_I, 1)).T
    return (L1(phi, tile_zeta, xi) - L2(phi, tile_zeta, xi)).sum()*V_I/N_I
    
def eval_L5(phi, zeta):
    xi, alpha, sin_theta = get_xi_alpha()
    N_I = xi.shape[1]
    #print N_I
    tile_zeta = np.tile(zeta, (N_I, 1)).T
    B = np.sum(alpha*(xi - tile_zeta), axis=0)
    zeta1, xi1 = tile_zeta + alpha*B, xi - alpha*B
    mask = (np.sum(np.square(zeta1/cut), axis=0) < 1) * (np.sum(np.square(xi1/cut), axis=0) < 1)
    #N_I = xi[:,mask].shape[1]
    #print N_I
    ci = phi(zeta1[:,mask]) + phi(xi1[:,mask]) - phi(zeta) - phi(xi[:,mask])
    return a3*(ci*np.exp(-sqr(xi[:,mask]))*np.abs(B[mask])*sin_theta[mask]).sum()*V_I/N_I
    #return a3*(ci*np.exp(-sqr(xi))*np.abs(B)*sin_theta).sum()*V_I/N_I

def eval_J(phi, psi, zeta):
    xi, alpha, sin_theta = get_xi_alpha()
    N_I = xi.shape[1]
    tile_zeta = np.tile(zeta, (N_I, 1)).T
    V = xi - tile_zeta
    B = np.sum(alpha*V, axis=0)
    zeta1, xi1 = tile_zeta + alpha*B, xi - alpha*B
    ci = phi(zeta1)*psi(xi1) + phi(xi1)*psi(zeta1) - phi(zeta)*psi(xi) - phi(xi)*psi(zeta)
    return .5*a3*(ci*np.exp(-sqr(xi))*np.abs(B)*sin_theta).sum()*V_I/N_I

def solve2(phi, Phi, get_zeta, f, f_old, F_old, alpha=1., pre_phi=lambda x: x[0]**0, beta=1):
    f_new = np.copy(f)
    F = np.copy(F_old)
    alpha_ = X**0
    for i in xrange(1,len(X)):
        zeta = np.transpose(get_zeta(X[i]))
        #print >> sys.stderr, "%.1e %+.5f" % (phi(zeta), delta)
        C = 1 + (pre_phi(zeta))**(-beta)
        #print >> sys.stderr, "pre: ", pre_phi(zeta)
        alpha_[i] = alpha / C
        if X[i] < 5:
            F[i] = (eval_L_plus_nu(phi, zeta) - Phi(zeta))/nu(mag(zeta))/pre_phi(zeta)
        else:
            F[i] = f[i] + (eval_L5(phi, zeta) - Phi(zeta)) / pre_phi(zeta)
    zero_corrector2(F, True)
    tau = (F - F_old) / (f - f_old) - 1
    print >> sys.stdout, (F-f)/f #,  alpha_ * np.sign(tau) * F
    #print >> sys.stderr, np.sign(tau)
    alpha_ *= -np.sign(tau) #*np.sign(F_old-f_old)
    f_new = alpha_ * F + (1.-alpha_) * f
    '''
    delta = eval_L3(phi, zeta) - Phi(zeta)
    print >> sys.stderr, mag(zeta), eval_L3(phi, zeta), Phi(zeta)
    f[i] += (delta)/(beta + sqr(zeta)**gamma) * alpha
    f_new[i] = alpha_ * (eval_L_plus_nu(phi, zeta) - Phi(zeta))/nu(mag(zeta))/pre_phi(zeta) + (1.-alpha_) * f[i]
    '''
    return f_new, f, F

def solve(pre_phi, Phi, get_zeta, f0, result, corrector = lambda f: None, alpha=1., beta=0.):
    f, F = f0(X), 0*X
    f_old = 0.99*f
    for n in xrange(int(sys.argv[2])):
        print n, f
        print >> sys.stderr, n, result(f)
        phi = lambda x: pre_phi(x)*interp1d(X, f)(mag(x))
        f, f_old, F = solve2(phi, Phi, get_zeta, f, f_old, F, alpha, pre_phi, beta)
        corrector(f)
    return f

def symmetrize(x, y):
    return np.append(-x, x), np.append(y, y)

def zero_corrector(f):
    XX, ff = symmetrize(X[1:-1], f[1:-1])
    f[0] = np.poly1d(np.polyfit(XX, ff, 32))(0)

def zero_corrector2(f, quiet=False, kind='quadratic'):
    func = {
        'linear': lambda x, a, b, c: a + b*x + c*x**2,
        'quadratic': lambda x, a, b: a + b*x**2,
        'cubic': lambda x, a, b, c: a + b*x**2 + c*x**3
    }[kind]
    try:
        f0new = curve_fit(func, X[1:4], f[1:4])[0][0]
    except RuntimeError:
        f0new = f[1]
    if not quiet:
        print >> sys.stderr, "  zero_corr2:", f0new - f[0]
    f[0] = f0new

def sign_corrector(f, sgn=1):
    f *= (np.sign(f*sgn)+1)/2

def eval_allL(phi, get_zeta, eval_L):
    f = 0*X
    for i in xrange(len(X)):
        zeta = np.transpose(get_zeta(X[i]))
        f[i] = eval_L(phi, zeta)
    return f

def eval_allJ(phi, psi, get_zeta):
    f = 0*X
    for i in xrange(len(X)):
        zeta = np.transpose(get_zeta(X[i]))
        f[i] = eval_J(phi, psi, zeta)
    return f

def eval_func(phi, get_zeta, begin=0):
    f = 0*X
    for i in xrange(begin,len(X)):
        zeta = np.transpose(get_zeta(X[i]))
        f[i] = phi(zeta)
    if begin > 0:
        zero_corrector2(f, True)
    return f

def average(func):
    total = 0*X
    N = int(sys.argv[2])
    for n in xrange(N):
        f = func()
        print n, f
        print >> sys.stderr, n, np.sum(f)
        total += f
    return total/N

def calc_I(power, f, beg=0):
    return 8./15/np.sqrt(np.pi) * np.trapz((X**power * np.exp(-X*X))[beg:]*f[beg:], X[beg:])

eta, A_, B_, D1_, D2_, F_ = np.loadtxt("../tables/ci-functions.txt").T

def calc_gamma1():
    B_0 = lambda x: x[0]*x[1]
    f0 = np.poly1d(np.polyfit(eta, B_, 16))
    #f1 = np.poly1d(np.polyfit(eta, B_+0.01, 16))
    eta_, f_, f1_, diff_ = np.loadtxt("./gamma1.txt").T
    f1 = np.poly1d(np.polyfit(eta_, f1_, 16))
    result = lambda f: calc_I(6, f)
    f = solve(B_0, lambda x: -2*B_0(x), zeta2d, f1, result, zero_corrector2, 0.5, 0.1)
    print >> sys.stderr, result(f), result(f0(X))
    py.plot(X, f0(X), 'r', X, f, 'g', X, f1(X), 'g--')
    show()
    np.savetxt(sys.stderr, np.transpose((X, f0(X), f, (f-f0(X))/f0(X))), fmt='%1.5e')

def calc_gamma2():
    A_0 = lambda x: x[2]
    def corrector(f):
        k = 3./8*np.sqrt(np.pi)
        C = np.trapz(X**4 * f * np.exp(-X*X), X)
        print >> sys.stderr, "  corr:", C
        f -= C / k
    f0 = np.poly1d(np.polyfit(eta, A_, 16))
    f1 = np.poly1d(np.polyfit(eta, A_, 16))
    result = lambda f: 2*calc_I(6, f)
    f = solve(A_0, lambda x: -A_0(x)*(sqr(x)-2.5), zeta1d, f1, result, corrector)
    print >> sys.stderr, result(f), result(f0(X))
    py.plot(X, f0(X), X, f)
    show()

def calc_gamma3():
    F_0 = lambda x: x[0]*x[1]
    f0 = np.poly1d(np.polyfit(eta, F_, 16))
    f1 = np.poly1d(np.polyfit(eta, F_+0.2, 16))
    A = np.poly1d(np.polyfit(eta, A_, 16))
    result = lambda f: -2*calc_I(6, f)
    f = solve(F_0, lambda x: F_0(x)*A(mag(x)), zeta2d, f1, result, zero_corrector)
    print >> sys.stderr, result(f), result(f0(X))
    py.plot(X, f0(X), X, f)
    show()

def calc_gamma9():
    K_0 = lambda x: x[0]*x[1]
    eta_, K_ = np.loadtxt("./gamma9.txt").T
    f0 = np.poly1d(np.polyfit(eta_, K_, 16))
    B = np.poly1d(np.polyfit(eta, B_, 16))
    result = lambda f: -calc_I(6, f)
    f = solve(K_0, lambda x: K_0(x)*B(mag(x)), zeta2d, f0, result, zero_corrector2, 0.25, 0.25)
    print >> sys.stderr, result(f), result(f0(X))
    py.plot(X, f0(X), '--', X, f)
    show()
    print >> sys.stderr, "--------------------------------------"
    np.savetxt(sys.stderr, np.transpose((X, f)), fmt='%1.5e')

def calc_gamma8a():
    m = np.arange(1, len(X))
    B = np.poly1d(np.polyfit(eta, B_, 16))
    derB = np.poly1d(np.polyfit(eta, B_, 16)).deriv()
    derB_ = np.vectorize(lambda x: -1.10 if x==0 else derB(x)/x)
    IB = lambda x: 2*B(x) - derB_(x)
    pre1 = lambda x: x[0]*x[1]
    pre2 = lambda x: (3*x[2]**2 - x[0]**2) * pre1(x)
    j_ = lambda x: 3*Q2_(mag(x)) + x[0]**2*Q3_(mag(x))
    phi1 = lambda x: pre1(x) * j_(x)
    phi2 = lambda x: pre2(x) * Q3_(mag(x))
    Phi1 = lambda x: -pre1(x) * IB(mag(x)) * x[0]**2
    Phi2 = lambda x: -pre2(x) * IB(mag(x))
    
    eta_, Q2_, Q3_ = np.loadtxt("./gamma8a.txt").T
    zero_corrector2(Q3_)
    kind = 'linear'
    #Q2_, Q3_ = interp1d(eta_, Q2_, kind=kind), interp1d(eta_, Q3_, kind=kind)
    Q2_, Q3_ = np.poly1d(np.polyfit(eta_, Q2_, 16)), np.poly1d(np.polyfit(eta_, Q3_, 16))
    Q2, Q3 = Q2_(X), Q3_(X)
    j = eval_func(j_, zeta2d, begin=1)
    j_old, Q3_old = 0.99*j, 0.99*Q3
    F_j, F_Q3 = 0*X, 0*X
    j_0, Q2_0, Q3_0 = np.copy(j), np.copy(Q2), np.copy(Q3)
    beta = 0.2
    for n in xrange(int(sys.argv[2])):
        Q2_, Q3_ = interp1d(X, Q2, kind=kind), interp1d(X, Q3, kind=kind)
        #Q3, Q3_old, F_Q3 = solve2(phi2, Phi2, zeta3d, Q3, Q3_old, F_Q3, 0.5, pre2, beta)
        j, j_old, F_j = solve2(phi1, Phi1, zeta2d, j, j_old, F_j, 0.2, pre1, beta)
        Q2[m] = (j[m] - .5*(X**2*Q3)[m]) / 3
        for f in [ j, Q2, Q3 ]:
            zero_corrector2(f)
        for f in [ j, Q2 ]:
            sign_corrector(f,1)
        result2 = lambda f: calc_I(6, f)
        result3 = lambda f: calc_I(8, f)/7
        print >> sys.stderr, n, result2(Q2), result3(Q3)
        print n
        #np.savetxt(sys.stdout, np.transpose((X, Q2, Q3)), fmt='%1.5e')
    
    #py.plot(X, j_0, 'r--', X, Q2_0, 'g--', X, Q3_0, 'b--', X, j, 'r', X, Q2, 'g', X, Q3, 'b')
    py.plot(X, Q2_0, 'g--', X, Q3_0, 'b--', X, Q2, 'g', X, Q3, 'b')
    show()
    print >> sys.stderr, "--------------------------------------"
    np.savetxt(sys.stderr, np.transpose((X, Q2, Q3)), fmt='%1.5e')

def calc_gamma8b():
    def correct_quattro(Y, x1, x2):
        fit = (X >= x1) * (X <= x2)
        corr = (X < x1)
        func = lambda x, a, b, c: (a + b*x**2 + c*x**3)*x**4
        Y[corr] = func(X[corr], *curve_fit(func, X[fit], Y[fit])[0])

    def interp_quattro(Y):
        temp = interp1d(np.append(0, X[1:]), np.append(0, Y[1:]/X[1:]**4))
        return lambda x: temp(x)*x**4
        
    m = np.arange(1, len(X))
    pre1 = lambda x: x[0]*x[1]
    pre2 = lambda x: 6*(x[0]*x[1])**2 - x[0]**4 - x[1]**4
    pre3 = lambda x: (3*x[2]**2 - x[0]**2) * pre1(x)
    pre4 = lambda x: pre2(x)-6*pre3(x)
    j1_ = lambda x: Q21_(mag(x)) + x[2]**2 * Q30_(mag(x))
    j2_ = lambda x: Q22_(mag(x)) + x[2]**2 * Q30_(mag(x))
    phi1 = lambda x: pre1(x) * j1_(x)
    phi2 = lambda x: pre1(x) * j2_(x)
    phi3 = lambda x: pre2(x) * Q30_(mag(x))
    phi4 = lambda x: pre3(x) * Q30_(mag(x))
    eta_, J1_, J2_, J3_, J4_, J5_, J6_ = np.loadtxt("./gamma8ih.txt").T
    J5, J6 = interp1d(eta_, J5_), interp1d(eta_, J6_)
    J123_ = J2_ + 2*J3_ - J1_
    J456_ = J5_ + 2*J6_ - J4_
    correct_quattro(J123_, .3, .6)
    correct_quattro(J456_, .2, .6)
    J123, J456 = interp_quattro(J123_), interp_quattro(J456_)
    Phi1 = lambda x: J5(mag(x))
    Phi2 = lambda x: J6(mag(x))
    Phi3 = lambda x: 2*J123(mag(x))
    Phi4 = lambda x: J456(mag(x))
    kind = 'linear'

    eta_, Q21_, Q22_, Q30_ = np.loadtxt("./gamma8b.txt").T
    #Q21_, Q22_, Q30_ = interp1d(eta_, Q21_, kind=kind), interp1d(eta_, Q22_, kind=kind), interp1d(eta_, Q30_, kind=kind)
    Q21_, Q22_, Q30_ = np.poly1d(np.polyfit(eta_, Q21_, 16)), np.poly1d(np.polyfit(eta_, Q22_, 16)), np.poly1d(np.polyfit(eta_, Q30_, 16))
    Q21, Q22, Q30 = Q21_(X), Q22_(X), Q30_(X)
    #Q21, Q22, Q3 = 0*X, 0*X, 0*X      ### COMMENT THIS
    j1, j2 = eval_func(j1_, zeta3d), eval_func(j2_, zeta3d)
    j1_old, j2_old, Q30_old = 0.99*j1, 0.99*j2, 0.99*Q30
    F_j1, F_j2, F_Q30 = 0*X, 0*X, 0 *X
    j1_0, j2_0, Q21_0, Q22_0, Q30_0 = np.copy(j1), np.copy(j2), np.copy(Q21), np.copy(Q22), np.copy(Q30)
    beta = 0.0
    for n in xrange(int(sys.argv[2])):
        Q21_, Q22_, Q30_ = interp1d(X, Q21, kind=kind), interp1d(X, Q22, kind=kind), interp1d(X, Q30, kind=kind)
        #Q30, Q30_old, F_Q30 = solve2(phi3, Phi3, zeta2d, Q30, Q30_old, F_Q30, 0.5, pre2, beta)
        Q30, Q30_old, F_Q30 = solve2(phi4, Phi4, zeta3d, Q30, Q30_old, F_Q30, 0.5, pre3, beta)
        #j1, j1_old, F_j1 = solve2(phi1, Phi1, zeta3d, j1, j1_old, F_j1, 0.5, pre1, beta)
        #j2, j2_old, F_j2 = solve2(phi2, Phi2, zeta3d, j2, j2_old, F_j2, 0.5, pre1, beta)
        Q21[m] = (j1[m] - (X**2*Q30)[m]/3)
        Q22[m] = (j2[m] - (X**2*Q30)[m]/3)
        for f in [ j1, j2, Q21, Q22 ]:
            zero_corrector2(f)
        for f in [ Q30 ]:
            zero_corrector2(f) #, kind='cubic')
        result2 = lambda f: calc_I(6, f)
        result3 = lambda f: calc_I(8, f, 1)/7
        print >> sys.stderr, n, result2(Q21), result2(Q22), result3(Q30)
        print n
        #np.savetxt(sys.stdout, np.transpose((X, Q22, Q30)), fmt='%1.5e')
    
    #py.plot(X, Q21_0, 'r--', X, Q22_0, 'g--', X, Q30_0, 'b--', X, Q21, 'r', X, Q22, 'g', X, Q30, 'b')
    py.plot(X, Q30_0, 'b--', X, Q30, 'b')
    show()
    print >> sys.stderr, "--------------------------------------"
    np.savetxt(sys.stderr, np.transpose((X, Q21, Q22, Q30)), fmt='%1.5e')

def calc_gamma8ih():
    B = np.poly1d(np.polyfit(eta, B_, 16))
    phi1 = lambda x: x[0]**2 * B(mag(x))
    phi2 = lambda x: x[1]**2 * B(mag(x))
    phi3 = lambda x: x[0]*x[1] * B(mag(x))
    phi4 = lambda x: x[0]*x[2] * B(mag(x))
    phi5 = lambda x: x[1]*x[2] * B(mag(x))
    eta_, J1_, J2_, J3_, J4_, J5_, J6_ = np.loadtxt("./gamma8ih.txt").T
    
    J1 = average(lambda: eval_allJ(phi1, phi1, zeta2d))
    J2 = average(lambda: eval_allJ(phi1, phi2, zeta2d))
    J3 = average(lambda: eval_allJ(phi3, phi3, zeta2d))
    J4 = average(lambda: eval_allJ(phi1, phi3, zeta3d))
    J5 = average(lambda: eval_allJ(phi1, phi5, zeta3d))
    J6 = average(lambda: eval_allJ(phi3, phi4, zeta3d))
    #J4, J5, J6 = J4_, J5_, J6_     ### COMMENT THIS
    J4[0] = J5[0] = J6[0] = 0
    py.plot(X, J1, 'r', eta_, J1_, 'r--', X, J2, 'g', eta_, J2_, 'g--', X, J3, 'b', eta_, J3_, 'b--')
    show()
    py.plot(X, J4, 'r', eta_, J4_, 'r--', X, J5, 'g', eta_, J5_, 'g--', X, J6, 'b', eta_, J6_, 'b--')
    show()
    np.savetxt(sys.stderr, np.transpose((X, J1, J2, J3, J4, J5, J6)), fmt='%1.5e')

def corrector_T(T1, T2, beta=1):
    k = 3./8*np.sqrt(np.pi)
    C = np.trapz(X**4 * (5*T1 + X**2*T2) * np.exp(-X*X), X)
    print >> sys.stderr, "  corr:", C
    T1 -= beta * C / k / 5

def calc_gamma10(letter):
    m = np.arange(1, len(X))
    A = np.poly1d(np.polyfit(eta, A_, 16))
    B = np.poly1d(np.polyfit(eta, B_, 16))
    derA = np.poly1d(np.polyfit(eta, A_, 16)).deriv()
    derA_ = np.vectorize(lambda x: 6.87 if x==0 else derA(x)/x)
    IA = lambda x: 2*A(x) - derA_(x)
    derB = np.poly1d(np.polyfit(eta, B_, 16)).deriv()
    IB = lambda x: (x**2-3)*B(x) - .5*x*derB(x)

    pre1 = lambda x: x[0]
    pre2 = lambda x: x[0]*x[1]*x[2]
    j_ = lambda x: 3*T1_(mag(x)) + x[0]**2*T2_(mag(x))
    phi1 = lambda x: pre1(x) * j_(x)
    phi2 = lambda x: pre2(x) * T2_(mag(x))
    if letter == 'a':
        Phi1 = lambda x: -pre1(x) * IA(mag(x)) * x[0]**2
        Phi2 = lambda x: -pre2(x) * IA(mag(x))
        eta_, T1_, T2_ = np.loadtxt("./gamma10a.txt").T
    else:
        Phi1 = lambda x: -pre1(x) * ( IB(mag(x)) * x[0]**2 - 1.5*gamma_1 )
        Phi2 = lambda x: -pre2(x) * IB(mag(x))
        eta_, T1_, T2_ = np.loadtxt("./gamma10b.txt").T
    
    #T1_ = 0*eta_ ### COMMENT THIS
    kind = 'linear'
    #T1_, T2_ = interp1d(eta_, T1_, kind=kind), interp1d(eta_, T2_, kind=kind)
    T1_, T2_ = np.poly1d(np.polyfit(eta_, T1_, 16)), np.poly1d(np.polyfit(eta_, T2_, 16))
    T1, T2 = T1_(X), T2_(X)
    j = eval_func(j_, zeta1d)
    j_old, T2_old = 0.99*j, 0.99*T2
    F_j, F_T2 = 0*X, 0*X
    j_0, T1_0, T2_0 = np.copy(j), np.copy(T1), np.copy(T2)
    beta = 0.0
    for n in xrange(int(sys.argv[2])):
        T1_, T2_ = interp1d(X, T1, kind=kind), interp1d(X, T2, kind=kind)
        #T1_, T2_ = np.poly1d(np.polyfit(X, T1, 16)), np.poly1d(np.polyfit(X, T2, 16))
        j = eval_func(j_, zeta1d)
        T2, T2_old, F_T2 = solve2(phi2, Phi2, zeta3d, T2, T2_old, F_T2, 0.5, pre2, beta)
        j, j_old, F_j = solve2(phi1, Phi1, zeta1d, j, j_old, F_j, 0.5, pre1, beta)
        T1[m] = (j[m] - (X**2*T2)[m]) / 3
        for f in [ j, T1, T2 ]:
            zero_corrector2(f)
        #for f in [ j, Q2 ]:
        #    sign_corrector(f,1)
        corrector_T(T1, T2)
        print >> sys.stderr, n, .625*calc_I(6, T1), .125*calc_I(8, T2)
        print n
    
    py.plot(X, T1_0, 'g--', X, T2_0, 'b--', X, T1, 'g', X, T2, 'b')
    show()
    print >> sys.stderr, "--------------------------------------"
    np.savetxt(sys.stderr, np.transpose((X, T1, T2)), fmt='%1.5e')

def calc_gamma10c():
    def interp_power(Y, n):
        temp = interp1d(np.append(0, X[1:]), np.append(0, Y[1:]/X[1:]**n))
        return lambda x: temp(x)*x**n
        
    m = np.arange(1, len(X))
    pre1 = lambda x: x[1]
    pre2 = lambda x: x[0]*x[1]*x[2]
    j_ = lambda x: T12_(mag(x)) + x[0]**2 * T2_(mag(x))
    phi1 = lambda x: pre1(x) * j_(x)
    phi2 = lambda x: pre2(x) * T2_(mag(x))
    eta_, J1_, J2_ = np.loadtxt("./gamma10ih.txt").T
    J1, J2 = interp_power(J1_, 1), interp_power(J2_, 3)
    Phi1 = lambda x: J1(mag(x))
    Phi2 = lambda x: J2(mag(x))
    kind = 'linear'

    eta_, T12_, T2_ = np.loadtxt("./gamma10c.txt").T
    #T12_, T2_ = interp1d(eta_, T12_, kind=kind), interp1d(eta_, T2_, kind=kind)
    T12_, T2_ = np.poly1d(np.polyfit(eta_, T12_, 16)), np.poly1d(np.polyfit(eta_, T2_, 16))
    T12, T2 = T12_(X), T2_(X)
    #T12, T2 = 0*X, 0*X      ### COMMENT THIS
    j = eval_func(j_, zeta2d)
    j_old, T2_old = 0.99*j, 0.99*T2
    F_j, F_T2 = 0*X, 0*X
    j_0, T12_0, T2_0 = np.copy(j), np.copy(T12), np.copy(T2)
    beta = 0.0
    for n in xrange(int(sys.argv[2])):
        T12_, T2_ = interp1d(X, T12, kind=kind), interp1d(X, T2, kind=kind)
        #T12_, T2_ = np.poly1d(np.polyfit(X, T12, 16)), np.poly1d(np.polyfit(X, T2, 16))
        j = eval_func(j_, zeta2d)
        #T2, T2_old, F_T2 = solve2(phi2, Phi2, zeta3d, T2, T2_old, F_T2, 0.5, pre2, beta)
        j, j_old, F_j = solve2(phi1, Phi1, zeta2d, j, j_old, F_j, 0.5, pre1, beta)
        T12[m] = (j[m] - (X**2*T2)[m]/2)
        for f in [ j, T12, T2 ]:
            zero_corrector2(f)
        corrector_T(T12, T2, 0.5)
        print >> sys.stderr, n, .625*calc_I(6, T12), .125*calc_I(8, T2)
        print n
    
    py.plot(X, T12_0, 'r--', X, T2_0, 'g--', X, T12, 'r', X, T2, 'g')
    show()
    print >> sys.stderr, "--------------------------------------"
    np.savetxt(sys.stderr, np.transpose((X, T12, T2)), fmt='%1.5e')

def calc_gamma10ih():
    A = np.poly1d(np.polyfit(eta, A_, 16))
    B = np.poly1d(np.polyfit(eta, B_, 16))
    phi1 = lambda x: x[0] * A(mag(x))
    phi2 = lambda x: x[0]*x[1] * B(mag(x))
    phi3 = lambda x: x[1]*x[2] * B(mag(x))
    eta_, J1_, J2_ = np.loadtxt("./gamma10ih.txt").T
    
    J1 = average(lambda: eval_allJ(phi1, phi2, zeta2d))
    J2 = average(lambda: eval_allJ(phi1, phi3, zeta3d))
    py.plot(X, J1, 'r', eta_, J1_, 'r--', X, J2, 'g', eta_, J2_, 'g--')
    show()
    np.savetxt(sys.stderr, np.transpose((X, J1, J2)), fmt='%1.5e')

def compare_func(phi, Phi, get_zeta, diff=True, pre=lambda x: X**0):
    L3 = eval_allL(phi, get_zeta, eval_L3)
    L5 = eval_allL(phi, get_zeta, eval_L5)
    F = eval_func(Phi, get_zeta)
    py.axhline(lw=.5, c='k', ls=':')
    if diff:
        py.plot(X, (L3-F)/pre(X), 'r', X, (L5-F)/pre(X), 'g')
    else:
        py.plot(X, L3/pre(X), 'r', X, L5/pre(X), 'g', X, F/pre(X), 'b--')
    show()

def compare_L3_L5():
    A = np.poly1d(np.polyfit(eta, A_, 16))
    B = np.poly1d(np.polyfit(eta, B_, 16))
    #one = lambda x: 0*x[0]+1
    #J1 = eval_allJ(phi, one, zeta1d)
    #J2 = eval_allJ(one, phi, zeta1d)
    '''
    compare_func(
        lambda x: x[0] * A(mag(x)),
        lambda x: -x[0] * (sqr(x)-2.5),
        zeta1d)
    compare_func(
        lambda x: x[0] * x[1] * B(mag(x)),
        lambda x: -2 * x[0] * x[1],
        zeta2d, False)
    '''
    '''
    # gamma9
    eta_, B4_ = np.loadtxt("./gamma9.txt").T
    B4 = interp1d(eta_, B4_)
    compare_func(
        lambda x: x[0] * x[1] * B4(mag(x)),
        lambda x: x[0] * x[1] * B(mag(x)),
        zeta2d)
    '''
    '''
    # gamma8a
    eta_, Q2_, Q3_ = np.loadtxt("./gamma8a.txt").T
    Q2_, Q3_ = interp1d(eta_, Q2_), interp1d(eta_, Q3_)
    derB = np.poly1d(np.polyfit(eta, B_, 16)).deriv()
    derB_ = np.vectorize(lambda x: -1.10 if x==0 else derB(x)/x)
    IB = lambda x: 2*B(x) - derB_(x)
    corr = lambda x: 2-np.log(x)
    compare_func(
        lambda x: x[0]*x[1]*(3*Q2_(mag(x)) + x[0]**2*Q3_(mag(x))*corr(mag(x))),
        lambda x: -x[0]**3*x[1]*IB(mag(x)),
        zeta2d, False,
        lambda x: x**3)
    compare_func(
        lambda x: x[1]*x[2]*(Q2_(mag(x)) + x[0]**2*Q3_(mag(x))*corr(mag(x))),
        lambda x: -x[0]**2*x[1]*x[2]*IB(mag(x)),
        zeta3d, False,
        lambda x: x**3)
    '''
    '''
    # gamma8b
    pre1, pre2 = lambda x: x[0]*x[1], lambda x: x[1]*x[2]
    j1_ = lambda x: Q21_(mag(x)) + 2*Q22_(mag(x)) + x[0]**2*Q30_(mag(x))
    j2_ = lambda x: Q21_(mag(x)) + x[0]**2 * Q30_(mag(x))
    j3_ = lambda x: Q22_(mag(x)) + x[0]**2 * Q30_(mag(x))
    phi1 = lambda x: pre1(x) * j1_(x)
    phi2 = lambda x: pre2(x) * j2_(x) 
    phi3 = lambda x: pre2(x) * j3_(x)
    eta_, J1_, J2_, J3_, J4_, J5_, J6_ = np.loadtxt("./gamma8ih.txt").T
    J4_, J5_, J6_ = interp1d(eta_, J4_), interp1d(eta_, J5_), interp1d(eta_, J6_)
    eta_, Q21_, Q22_, Q30_ = np.loadtxt("./gamma8b.txt").T
    Q21_, Q22_, Q30_ = interp1d(eta_, Q21_), interp1d(eta_, Q22_), interp1d(eta_, Q30_)
    compare_func(phi1, lambda x: J4_(mag(x)), zeta2d, False)
    compare_func(phi2, lambda x: J5_(mag(x)), zeta3d, False)
    compare_func(phi3, lambda x: J6_(mag(x)), zeta3d, False)
    '''

### ----------------------- ###

{
    '1': calc_gamma1,
    '2': calc_gamma2,
    '3': calc_gamma3,
    '8a': calc_gamma8a,
    '8b': calc_gamma8b,
    '8ih': calc_gamma8ih,
    '9': calc_gamma9,
    '10a': partial(calc_gamma10, 'a'),
    '10b': partial(calc_gamma10, 'b'),
    '10c': calc_gamma10c,
    '10ih': calc_gamma10ih,
    'compare': compare_L3_L5
}[sys.argv[1]]()

