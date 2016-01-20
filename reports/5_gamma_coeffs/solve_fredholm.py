#!/usr/bin/env python

import sys
import os
import numpy as np
#import pylab as py
from scipy.special import erf
from scipy.misc import comb
from datetime import datetime
from numpy.linalg import solve
from functools import partial
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

xi_min, xi_cut = 0., 6.4
_, problem, dirname = sys.argv
stride = int(dirname)
N = 100*stride
if not os.path.exists(dirname):
    print "Create directory %s, N = %d" % (dirname, N)
    os.mkdir(dirname)
h = (xi_cut-xi_min)/N
#xi_points = np.linspace(xi_min, xi_cut, N+1)
#xi = 0.5*(xi_points[0:N] + xi_points[1:N+1])
xi = np.linspace(xi_min, xi_cut, N+1)

k = lambda i: i.reshape(-1,1,1)
R_s_ = lambda x, y: x**2 + y**2
R_d_ = lambda x, y: np.fabs(x**2 - y**2)
r_s_ = lambda x, y: x + y
r_d_ = lambda x, y: np.fabs(x - y)
Erf = lambda x: np.sqrt(np.pi)/2*erf(x)

def calc_C(n, x, y, C):
    R_s, R_d = R_s_(x,y), R_d_(x,y)
    r_s, r_d = r_s_(x,y), r_d_(x,y)
    return {
        0: lambda: 2*np.exp((R_s-R_d)/2)*Erf((r_s-r_d)/2),
        1: lambda: (2+R_d)*C[0] - 2*(r_s-r_d)
    }.get(n, lambda: 2*(2*n-1)*C[n-1] + R_d**2*C[n-2] - 2*(r_s**(2*n-1)-r_d**(2*n-1)))()

def G(n, x, y):
    R_s = R_s_(x,y)
    C = np.zeros((n+1, x.size, y.size))
    for i in xrange(n+1):
        C[i] = calc_C(i, x, y, C)
    a = lambda i: comb(n,k(i)) * R_s**(n-k(i)) * (-1)**k(i) * C[i]
    A = np.sum(np.fromfunction(a, (n+1,), dtype=int), axis=0)
    return A / (2**n*(x*y)**(n+1))

def J(n, x, y):
    R_s, r_s, r_d = R_s_(x,y), r_s_(x,y), r_d_(x,y)
    a = lambda i: comb(n,k(i)) * R_s**(n-k(i)) * (-1)**k(i) / (2*k(i)+3) * (r_s**(2*k(i)+3) - r_d**(2*k(i)+3))
    A = np.sum(np.fromfunction(a, (n+1,), dtype=int), axis=0)
    return A / (2**n*(x*y)**(n+1))

def interp(X, Y, n=0, deg=16, spline=False):
    if n:
        parabola = lambda x, a, b: a + b*x**2
        Y[1:] = Y[1:]/X[1:]**n
        Y[0] = curve_fit(parabola, X[1:4], Y[1:4])[0][0]
    if spline:
        return interp1d(X, Y, kind='cubic')
    else:
        return np.poly1d(np.polyfit(X, Y, deg))

def read(filename, pos=0, ipow=0, deg=16, value=None):
    try:
        _, f = np.loadtxt(os.path.join(dirname, filename + '.txt')).T
        print "read %s from %s" % (filename, dirname)
        return lambda x: f
    except IOError:
        try:
            data = np.loadtxt('tables/%s.txt' % filename).T
        except IOError:
            data = np.loadtxt('data/%s.txt' % filename).T
    pos = 1 if len(data) == 2 else pos
    spline = 'ihe' in filename
    if value:
        return interp(data[0], value(data), n=ipow, deg=deg, spline=spline)
    if pos:
        return interp(data[0], data[pos], n=ipow, deg=deg, spline=spline)
    else:
        return [ interp(data[0], data[i], n=ipow, deg=deg, spline=spline) for i in xrange(1, len(data)) ]

nu = lambda x: np.where(x > 0, np.exp(-x**2) + (2*x+1/x) * Erf(x), 2.)/2**1.5
#calc_I = lambda n, xi, f: 8./15/np.sqrt(np.pi)*np.sum(f*xi**n*np.exp(-xi**2))*h
calc_I = lambda n, xi, f: 8./15/np.sqrt(np.pi)*np.trapz(f*xi**n*np.exp(-xi**2), xi)
zero_corrector = lambda xi, f: 0

def plot_kernel(xi, Z):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    Z[np.diag_indices(len(Z))] += nu(xi)
    '''
    X, Y = np.meshgrid(xi, xi, indexing='ij')
    ax.plot_wireframe(X, Y, Z, rstride=stride, cstride=stride)
    '''
    x_max, y_max = 0.05, 0.05
    X, Y = np.meshgrid(xi[xi<x_max], xi[xi<y_max], indexing='ij')
    XX, YY = np.meshgrid(xi, xi, indexing='ij')
    m = (XX < x_max) & (YY < y_max)
    ax.plot_wireframe(X, Y, Z[m].reshape(X.shape))
    #'''
    ax.set_xlabel('x(zeta)')
    ax.set_ylabel('y(xi)')
    ax.set_zlabel('K(x,y)')
    py.show()

def rule_trap():
    rule = np.ones(xi.size)
    rule[0] = rule[-1] = 0.5
    return rule

### just use default values
def create(K, phi, func, value, corrector=zero_corrector, ipow=6, coeff=1.):
    return K, phi, func, value, corrector, ipow, coeff

### int_0^\infty K(x,y)f(y)dy - \nu(x)f(x) = phi(x)
def solve_ieqn(problem):
    K, phi, expected_func, expected_value, corrector, ipow, coeff = problems[problem]
    Z = np.zeros((xi.size, xi.size))
    print '[%s]' % datetime.now().time(), "Calculations are started"
    Z = K(xi.reshape(-1, 1), xi)*rule_trap()
    Z[np.diag_indices(xi.size)] -= nu(xi)
    print '[%s]' % datetime.now().time(), "Kernel has been calculated"
    f = solve(Z, phi(xi))
    print '[%s]' % datetime.now().time(), "Linear equations have been solved", np.allclose(np.dot(Z, f), phi(xi))
    if corrector(xi, 0*xi+1) != 0:
        C = corrector(xi, 0*xi+1) - corrector(xi, 0*xi)
        print 'correction =', corrector(xi, f)/C
        f -= corrector(xi, f)/C
    result = coeff*calc_I(ipow, xi, f)
    print "%g*I_%d = %.12f %.2e" % (coeff, ipow, result, np.fabs(result-expected_value))
    #py.plot(xi, f, 'r-', xi, expected_func(xi), 'g-')
    #py.show()
    np.savetxt(os.path.join(dirname, problem.upper() + '.txt'), np.transpose((xi, f)), fmt='%.10f')
    #np.savetxt(os.path.join(dirname, problem.upper()), np.transpose((xi[::stride], f[::stride])), fmt='%.10f')
    #plot_kernel(xi, Z)

pre = lambda y: np.exp(-y**2) * h / 2**1.5
preK = lambda n,x,y,F: np.where(y==0, 0, y**(2+n)/x**n * F)
asymK_1 = lambda x,y: 4*y/3*( (2 + y**2) + 11*x**2/5 + 8*x**4/7 + 16*x**6/45 )
asymK_2 = lambda x,y: 4*y/5*( (2 - y**2) + x**2*(-16*y**2/21 + 23./7) + x**4*(-16*y**2/63 + 40./21) + x**6*(-128*y**2/2079 + 64./99) )
asymK_3 = lambda x,y: 4*y/7*( (2 - 11*y**2/5) + x**2*(-16*y**2/9 + 13./3) + x**4*(-112*y**2/165 + 280./99) + x**6*(-128*y**2/715 + 448./429) )
asymK_4 = lambda x,y: 4*y/9*( (2 - 23*y**2/7 + 8*y**4/35) + x**2*(48*y**4/385 - 240*y**2/77 + 59./11)
    + x**4*(192*y**4/5005 - 192*y**2/143 + 560./143) + x**6*(128*y**4/15015 - 384*y**2/1001 + 224./143) )

asymK_3_ = lambda x,y: 4*y/7*( 2*(y/x)**7 + y**2*(13*(y/x)**7/3 - 11*(y/x)**5/5) + y**4*(280*(y/x)**7/99 - 16*(y/x)**5/9) + y**6*(448*(y/x)**7/429 - 112*(y/x)**5/165) )
asymK_4_ = lambda x,y: 4*y/9*( 2*(y/x)**9 + y**2*(59*(y/x)**9/11 - 23*(y/x)**7/7) + y**4*(560*(y/x)**9/143 - 240*(y/x)**7/77 + 8*(y/x)**5/35)
    + y**6*(224*(y/x)**9/143 - 192*(y/x)**7/143 + 48*(y/x)**5/385) )

K_1_ = lambda x,y: preK(1,x,y, 4*G(1,x,y)-2*J(1,x,y))
K_2_ = lambda x,y: preK(2,x,y, 6*G(2,x,y)-2*G(0,x,y)-3*J(2,x,y)+J(0,x,y))
K_3_ = lambda x,y: preK(3,x,y, 10*G(3,x,y)-6*G(1,x,y)-5*J(3,x,y)+3*J(1,x,y))
K_4_ = lambda x,y: preK(4,x,y, (70*G(4,x,y)-60*G(2,x,y)+6*G(0,x,y)-35*J(4,x,y)+30*J(2,x,y)-3*J(0,x,y))/4)

K_1 = lambda x,y: pre(y)*np.where((x<1e-3)&(y>=x), asymK_1(x,y), K_1_(x,y))
K_2 = lambda x,y: pre(y)*np.where((x<5e-3)&(y>=x), asymK_2(x,y), K_2_(x,y))
K_3 = lambda x,y: pre(y)*np.where((x<5e-2)&(y>=x), asymK_3(x,y), np.where((x<2e-2)&(y<x), asymK_3_(x,y), K_3_(x,y)))
K_4 = lambda x,y: pre(y)*np.where((x<1.2e-1)&(y>=x), asymK_4(x,y), np.where((x<8e-2)&(y<x), asymK_4_(x,y), K_4_(x,y)))

eta, A_, B_, D1_, D2_, F_ = np.loadtxt("./tables/ci-functions.txt").T
D1, D2, F = interp(eta, D1_), interp(eta, D2_), interp(eta, F_)
asym_A = np.poly1d([ 0.18009442, 0, -0.58584376, 0, 3.43757195, 0, -6.13706449 ])
asym_B = np.poly1d([ -0.05323961, 0, 0.15546985, 0, -0.55094231, 0, 3.51038226 ])
### read table data
A, B = read('A'), read('B')
derA_, derB_ = np.gradient(A(0), h, edge_order=2), np.gradient(B(0), h, edge_order=2)
derA = lambda x: np.where(x < 0.5, asym_A.deriv()(x), derA_)
derB = lambda x: np.where(x < 0.5, asym_B.deriv()(x), derB_)
derA_x = lambda x: np.where(x < 0.1, np.poly1d(asym_A.deriv().c[:-1])(x), derA_/xi)
derB_x = lambda x: np.where(x < 0.1, np.poly1d(asym_B.deriv().c[:-1])(x), derB_/xi)
'''
### interpolate data from literature
num_A, num_B = read('A', deg=40), read('B', deg=30)
asym_A = np.poly1d([ 0.18009442, 0, -0.58584376, 0, 3.43757195, 0, -6.13706449 ])
asym_B = np.poly1d([ -0.05323961, 0, 0.15546985, 0, -0.55094231, 0, 3.51038226 ])
A = lambda x: np.where(x < 0.25, asym_A(x), num_A(x))
B = lambda x: np.where(x < 0.25, asym_B(x), num_B(x))
derA = lambda x: np.where(x < 0.25, asym_A.deriv()(x), num_A.deriv()(x))
derB = lambda x: np.where(x < 0.25, asym_B.deriv()(x), num_B.deriv()(x))
derA_x = lambda x: np.where(x < 0.25, np.poly1d(asym_A.deriv().c[:-1])(x), num_A.deriv()(x)/x)
derB_x = lambda x: np.where(x < 0.25, np.poly1d(asym_B.deriv().c[:-1])(x), num_B.deriv()(x)/x)
'''

gamma_1  = 1.270042427
gamma_2  = 1.922284066
gamma_3  = 1.947906335
gamma_7a = -0.317776895
gamma_7b = 0.128572890
gamma_9  = 1.636075500

gamma_t1_1  = 0.984935377
gamma_t1_2  = 0.505576785
gamma_t2_1  = 0.507277037
gamma_t2_2  = 0.540997311
gamma_tt_12 = 0.037652844
gamma_tt_2  = 0.006842638

gamma_q_2   = 0.54475338
gamma_q_3   = 0.993545571
gamma_qq_22 = 0.121142330
gamma_qq_3  = -0.078789047

print 'gamma_7 =', -(gamma_7a + gamma_7b)
print 'gamma_8 =', (gamma_q_2 - gamma_qq_22) + (gamma_q_3 - gamma_qq_3)
print 'gamma_10 =', (gamma_t1_1 + gamma_t2_1 - 2*gamma_tt_12) + (gamma_t1_2 + gamma_t2_2 - 2*gamma_tt_2)

T1_1, T1_2 = read('gamma10a')
T2_1, T2_2 = read('gamma10b')
TT_12, TT_2 = read('gamma10c')
Q_2, Q_3 = read('gamma8a') 
QQ_21, QQ_22, QQ_3 = read('gamma8b')

T0_2 = read('T0_2')
T1_2 = read('T1_2')
T2_2 = read('T2_2')
TT_2 = read('TT_2')
Q_3 = read('Q_3')
QQ_3 = read('QQ_3')

IT_1 = lambda x: 2*A(x) - derA_x(x)
IT_2 = lambda x: (x**2-3)*B(x) - x/2*derB(x)
ITT_3 = read('gamma10ih', value=lambda J: J[3] - 2*J[4], ipow=3)
ITT_1 = read('gamma10ih', value=lambda J: J[3] + 3*J[4], ipow=1)

IQ = lambda x: 2*B(x) - derB_x(x)
# gamma8ih2: zzzz, xxxx, xzxz, xxzz, xyxy, xxyy
IQQ_4 = read('gamma8ih2', value=lambda J: 4*J[1] + 3*J[2] - 16*J[3] - 8*J[4] + 2*J[5] + J[6], ipow=4)
IQQ_2 = read('gamma8ih2', value=lambda J: J[1] - J[2] + 3*J[3] - 2*J[4] - 3*J[5] + 2*J[6], ipow=2)

problems = {
    'a': create( K_1, lambda x: -(x**2-2.5), lambda x: A(x), gamma_2, corrector=partial(calc_I,4), coeff=2. ),
    'b': create( K_2, lambda x: 0*x-2, B, gamma_1 ),
    'b1': create( K_2, lambda x: -A(x), lambda x: -F(x), gamma_3, coeff=2. ),
    'b2': create( K_2, lambda x: 2*(x**2-3)*A(x)-x*derA(x), read('gamma7a'), gamma_7a ),
    'b3': create( K_2, lambda x: 4*read('gamma7ih', 1, ipow=2)(x), read('gamma7b'), gamma_7b ),
    'b4': create( K_2, lambda x: -B(x), lambda x: -read('gamma9')(x), gamma_9 ),
    't0_2': create( K_3, lambda x: -B(x), D2, calc_I(8, xi, D2(xi)), ipow=8 ),
    't0_1': create( K_1, lambda x: gamma_1 - 0.2*(x**2*B(x) + np.dot(K_1(x.reshape(-1, 1),xi)*xi**2, T0_2(xi)) - x**2*nu(x)*T0_2(x)), D1, calc_I(6, xi, D1(xi) ),
        corrector=lambda xi,f: 5*calc_I(4,xi,f) + calc_I(6,xi,D2(xi)) ),
    't1_2': create( K_3, lambda x: -IT_1(x), T1_2, gamma_t1_2, ipow=8, coeff=.125 ),
    't1_1': create( K_1, lambda x: -0.2*(x**2*IT_1(x) + np.dot(K_1(x.reshape(-1, 1),xi)*xi**2, T1_2(xi)) - x**2*nu(x)*T1_2(x)), T1_1, gamma_t1_1,
        corrector=lambda xi,f: 5*calc_I(4,xi,f) + calc_I(6,xi,T1_2(xi)), coeff=.625 ),
    't2_2': create( K_3, lambda x: -IT_2(x), T2_2, gamma_t2_2, ipow=8, coeff=.125 ),
    't2_1': create( K_1, lambda x: gamma_1/2 - 0.2*(x**2*IT_2(x) + np.dot(K_1(x.reshape(-1, 1),xi)*xi**2, T2_2(xi)) - x**2*nu(x)*T2_2(x)), T2_1, gamma_t2_1,
        corrector=lambda xi,f: 5*calc_I(4,xi,f) + calc_I(6,xi,T2_2(xi)), coeff=.625 ),
    'tt2': create( K_3, lambda x: ITT_3(x), TT_2, gamma_tt_2, ipow=8, coeff=.125 ),
    'tt12': create( K_1, lambda x: 0.2*(ITT_1(x) - np.dot(K_1(x.reshape(-1, 1),xi)*xi**2, TT_2(xi)) + x**2*nu(x)*TT_2(x)), TT_12, gamma_tt_12,
        corrector=lambda xi,f: 5*calc_I(4,xi,f) + calc_I(6,xi,TT_2(xi)), coeff=.625 ),
    'q3': create( K_4, lambda x: -IQ(x), Q_3, gamma_q_3, ipow=8, coeff=1./7 ),
    'q2': create( K_2, lambda x: -(x**2*IQ(x) + np.dot(K_2(x.reshape(-1, 1),xi)*xi**2, Q_3(xi)) - x**2*nu(x)*Q_3(x))/7, Q_2, gamma_q_2 ),
    'qq3': create( K_4, lambda x: IQQ_4(x)/4, QQ_3, gamma_qq_3, ipow=8, coeff=1./7 ),
    'qq22': create( K_2, lambda x: (IQQ_2(x) - np.dot(K_2(x.reshape(-1, 1),xi)*xi**2, QQ_3(xi)) + x**2*nu(x)*QQ_3(x))/7, QQ_22, gamma_qq_22 ),
}

with np.errstate(divide='ignore', invalid='ignore'):
    solve_ieqn(problem.lower())
