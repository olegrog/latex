#!/usr/bin/env python
r"Script for evaluating the errors in the various moments calculation of the discrete distribution function"
r"Usage: ./discrete_error.py <U> <N_r> <N_y> <cut_r> <ratio>"
r"Example:  ./discrete_error.py .5 12 22 4.3 1.3"

import numpy as np
from functools import partial
import sys

Rho, Temp, Speed = 1., 1., np.array([0,0,0],dtype=float)
Qflow, Tau = np.array([0.1, 0.2, 0]), np.array([0.05, 0., 0.2])
temp_ratio = 1.1    # let be > 1
Speed[0] = float(sys.argv[1])
N_r = int(sys.argv[2])
N_y = int(sys.argv[3])
heat = 0.1*(2*Speed[0])**2
cut_r = float(sys.argv[4]) # * np.sqrt(Temp*(1+heat))
q = float(sys.argv[5])
cut_new = 0.1*int(10*cut_r*np.sqrt(Temp*(1+heat)))
cut_x, cut_r = max(cut_r + Speed[0], cut_new), cut_new
N_x, N_z = int(N_r*cut_x/cut_r), N_r - 2 - int(2*Speed[0])

Tauij = [[0,Tau[2],Tau[1]],[Tau[2],0,Tau[0]],[Tau[1],Tau[0],0]]
zeros, ones = np.full(3, 0), np.full(3, 1)
dimension = (2*N_x, 2*N_y, 2*N_z)

def hermite(cut, N):
    H_xi, H_h = np.polynomial.hermite.hermgauss(2*N)
    H_xi *= np.sqrt(Temp)
    H_h /= np.exp(-H_xi**2)
    if np.sum(H_h) < 2*cut:
        C = (2*cut) / np.sum(H_h)
        H_xi *= C
        H_h *= C
    return np.vectorize(lambda i: H_xi[i]), np.vectorize(lambda i: H_h[i])

# width of the first cell
h1 = lambda q: cut_r/N_y if q==1 else cut_r*(q-1)/(q**N_y-1)

def i2h(i, q):
    j = lambda i: abs(i-N_y+.5)+.5
    nonuniform = lambda i: h1(q)*q**(j(i)-1)
    uniform = lambda i: cut_x/N_x + i*0
    return {
        0: hermite(cut_r, N_z)[1],  # hermite grid
        -1: uniform                 # uniform grid
    }.get(q, nonuniform)(i)         # exp refinement

def i2xi(i, q):
    sgn = lambda i: np.sign(i-N_y+.5)
    j = lambda i: abs(i-N_y+.5)+.5
    j_x = lambda i: i-N_x+.5
    nonuniform = lambda i: sgn(i)*h1(q) * (j(i)-.5 if q==1 else (q**j(i)+q**(j(i)-1)-2)/(q-1)/2)
    uniform = lambda i: j_x(i)*cut_x/N_x
    return {
        0: hermite(cut_r, N_z)[0],  # hermite grid
        -1: uniform                 # uniform grid
    }.get(q, nonuniform)(i)         # exp refinement

dxi = lambda i,j,k: i2h(i,-1)*i2h(j,q)*i2h(k,0)
xi = lambda i,j,k,v=zeros: (i2xi(i,-1)-v[0], i2xi(j,q)-v[1], i2xi(k,0)-v[2])
xi_x = lambda i,j,k,v=zeros: i2xi(i,-1)-v[0]
xi_y = lambda i,j,k,v=zeros: i2xi(j,q)-v[1]
sqr_xi = lambda i,j,k,v=zeros,s=ones: (s[0]*(i2xi(i,-1)-v[0]))**2 + (s[1]*(i2xi(j,q)-v[1]))**2 + (s[2]*(i2xi(k,0)-v[2]))**2
Maxwell = lambda v,t,i,j,k: Rho/(np.pi*t)**1.5 * np.exp(-sqr_xi(i,j,k,v)/t)*dxi(i,j,k)

# Cutting prolate spheroid 
radii = np.array([cut_x, cut_r, cut_r])

ind = np.fromfunction(lambda i,j,k: sqr_xi(i,j,k,zeros,radii**-1) <= 1, dimension)
u = np.fromfunction(lambda i,j,k: xi(i,j,k), dimension)
Dxi = np.fromfunction(lambda i,j,k: dxi(i,j,k), dimension)

X_x = np.fromfunction(lambda i: i2xi(i,-1), (2*N_x,))
X_y = np.fromfunction(lambda i: i2xi(i,q), (2*N_y,))
X_z = np.fromfunction(lambda i: i2xi(i,0), (2*N_z,))
H_x = np.fromfunction(lambda i: i2h(i,-1), (2*N_x,))
H_y = np.fromfunction(lambda i: i2h(i,q), (2*N_y,))
H_z = np.fromfunction(lambda i: i2h(i,0), (2*N_z,))

print_grid = False
if print_grid:
    print("--- x:", X_x, H_x)
    print("--- y:", X_y, H_y)
    print("--- z:", X_z, H_z)
#print 0.5*(H_y[1:]+H_y[:-1])   # should coincide
def calc_sizes(N,q):
    h_min, h_max = i2h(N,q), i2h(0,q)
    return (h_min, h_max, h_max/h_min)
print("Cut_x = %g, cut_r = %g, ratio = %g, N_z = %d" % (cut_x, cut_r, q, N_z))
print("Cell size (x): min = %.4g, max = %.4g, ratio = %.3g" % calc_sizes(N_x, -1))
print("Cell size (y): min = %.4g, max = %.4g, ratio = %.3g" % calc_sizes(N_y, q))
print("Cell size (z): min = %.4g, max = %.4g, ratio = %.3g" % calc_sizes(N_z, 0))
total = lambda X, cut: np.sum(np.abs(X) <= cut)/2
print("Total cells: %d (%d, %d, %d):" % (np.sum(ind), total(X_x, cut_x), total(X_y, cut_r), total(X_z, cut_r)))

def splot(f):
    import pylab as py
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = py.figure()
    ax = p3.Axes3D(fig)
    Xi = np.fromfunction(lambda i: i2xi(i,-1), (2*N_x,))
    Yi = np.fromfunction(lambda i: i2xi(i,q), (2*N_y,))
    xv, yv = np.meshgrid(Xi, Yi, sparse=False, indexing='ij')
    f[np.invert(ind[:,:,N_z])] = np.NaN
    f[ind[:,:,N_z]] /= Dxi[ind[:,:,N_z]]
    #ax.plot_wireframe(xv, yv, np.log(f[:,:,N_z]), linewidth=0.25)
    ax.plot_wireframe(xv, yv, f[:,:,N_z], linewidth=0.25)
    py.savefig('tmp.pdf')
    py.show()

def err(theor, real):
    return np.vectorize(lambda x,y: np.NaN if x==0 else abs(x-y)/x)(theor, real)

def calc_macro(f):
    rho = np.sum(f[ind])
    speed, qflow, tau = [0,0,0], [0,0,0], [0,0,0]
    for i in range(3):
        speed[i] = np.sum((f*u[i])[ind])/rho
    c = np.fromfunction(lambda i,j,k: xi(i,j,k,speed), dimension)
    sqr_c = np.fromfunction(lambda i,j,k: sqr_xi(i,j,k,speed), dimension)
    temp = 2./3*np.sum((f*sqr_c)[ind])/rho
    for i in range(3):
        qflow[i] = np.sum((f*c[i]*sqr_c)[ind])
        tau[i] = np.sum(2*(f*c[(i+1)%3]*c[(i+2)%3])[ind])
    return rho, temp, speed, qflow, tau

def test_1(Temp, Speed):
    P = Rho * Temp
    f1 = 1./P/Temp * np.fromfunction(lambda i,j,k: np.sum(np.tensordot(Tauij, xi(i,j,k,Speed), (0,0)) * xi(i,j,k,Speed), 0), dimension)
    f2 = 4./5/P/Temp * np.fromfunction(lambda i,j,k: np.tensordot(Qflow, xi(i,j,k,Speed), 1)*(sqr_xi(i,j,k,Speed)/Temp-2.5), dimension)
    f = np.fromfunction(partial(Maxwell, Speed, Temp), dimension) * (1 + f1 + f2)
    rho, temp, speed, qflow, tau = calc_macro(f)
    #splot(f)
    print("\n-- Test #1: continual flows - Grad's 13-moment approximation (temp = %g, speed = %g)" % (Temp, Speed[0]))
    print("rho =", err(Rho, rho))
    print("temp =", err(Temp, temp))
    print("speed =", err(Speed, speed))
    print("qflow =", err(Qflow, qflow/rho))
    print("tau =", err(Tau, tau/rho))

def test_2():
    Temp1, Temp2 = Temp, Temp/temp_ratio
    double_temp = np.sqrt(Temp1*Temp2)
    ss_temp = np.sqrt(Temp1) + np.sqrt(Temp2)
    delta_temp = np.abs(Temp2 - Temp1)

    Rho1, Rho2 = 2*Rho*np.sqrt(Temp2)/ss_temp, 2*Rho*np.sqrt(Temp1)/ss_temp
    f = Rho1/(np.pi*Temp1)**1.5 * np.fromfunction(lambda i,j,k: np.exp(-sqr_xi(i,j,k,Speed)/(Temp1))*dxi(i,j,k), dimension)
    negative = np.fromfunction(lambda i,j,k: j<N_y, dimension)
    f[negative] = Rho2/(np.pi*Temp2)**1.5 * np.fromfunction(lambda i,j,k: np.exp(-sqr_xi(i,j,k,-Speed)/(Temp2))*dxi(i,j,k), dimension)[negative]

    rho, temp, speed, qflow, tau = calc_macro(f)

    Rho_ = Rho
    Temp_ = double_temp * (1 + 8./3*np.dot(Speed,Speed)/ss_temp**2)
    Speed_ = -Speed * delta_temp / ss_temp**2
    Qflow_ = [ 0, 2*Rho*double_temp*delta_temp/ss_temp/np.sqrt(np.pi) * (1 + 2*np.dot(Speed,Speed)/ss_temp**2), 0 ]
    Tau_ = [ 0, 0, 4*Rho*Speed[0]*double_temp/ss_temp/np.sqrt(np.pi) ]

    #splot(f)
    #splot(f*c[0]*c[1])

    print("\n-- Test #2: free molecular flows - sum of 2 half-Maxwellians")
    print("rho =", err(Rho_, rho))
    print("temp =", err(Temp_, temp))
    print("speed =", err(Speed_, speed))
    print("qflow =", err(Qflow_, qflow/rho))
    print("tau =", err(Tau_, tau/rho))

def test_3():
    kn = 1
    gamma_1 = 1.270042427
    eta_, A_, B_, D1_, D2_, F_ = np.loadtxt("../tables/ci-functions.txt").T
    B = np.poly1d(np.polyfit(eta_, B_, 16))
    phi = np.fromfunction(lambda i,j,k: 2*Speed[0]*xi_x(i,j,k)*(1 - xi_y(i,j,k)*B(np.sqrt(sqr_xi(i,j,k)))*kn), dimension)
    f = np.fromfunction(partial(Maxwell, zeros, Temp), dimension) * (1 + phi)
    rho, temp, speed, qflow, tau = calc_macro(f)

    print("\n-- Test #3: linear case - asymptotic solution for hard-sphere molecules")
    print("rho =", err(Rho, rho))
    print("speed =", err(Speed, speed))
    print("tau =", err((0, 0, -2*gamma_1*Speed[0]*kn), tau/rho))

test_1(Temp, Speed)
test_1(Temp*(1 + heat), np.sqrt(Speed)/40)
test_2()
test_3()
