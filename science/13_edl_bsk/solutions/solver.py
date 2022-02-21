#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy import sinh, cosh, sqrt, log

parser = argparse.ArgumentParser(
    description='Solver for analysis of electric double layer with surface roughness')

str2pair = lambda s: [float(item) for item in s.split(':')]

parser.add_argument('-N', type=int, default=100, help='number of points')
parser.add_argument('-P', '--Psi', type=float, default=10, help='potential of the electrode')
parser.add_argument('-a', '--alpha', type=float, default=1, help='(sigma^2/lambda_c/lambda_D)^2')
parser.add_argument('-e', '--epsilon', type=float, default=1, help='sigma/lambda_c')
parser.add_argument('-g', '--gamma', type=float, default=1, help='volume fraction of ions')
parser.add_argument('-s', '--figsize', type=str2pair, default='5:4', help='figure size')
parser.add_argument('--pdf', action='store_true', help='save a PDF file instead')
args = parser.parse_args()

class Style:
    dotted = { 'linestyle': ':', 'linewidth': 0.5, 'color': 'k' }
    dashed = { 'linestyle': '--', 'linewidth': 0.5, 'color': 'k' }
    thin = { 'linestyle': '-', 'linewidth': 0.5, 'color': 'k' }

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
        'figure.figsize': args.figsize
    }
    plt.rcParams.update(params)

rhoC = args.alpha/args.epsilon**4
rho = lambda psi: sinh(psi)/(1 + 2*args.gamma*sinh(psi/2)**2)
int_rho = lambda psi: log(args.gamma*cosh(psi) - args.gamma + 1)/args.gamma
f = lambda z, psi: [ psi[1], psi[2], psi[3], psi[2] - rhoC*rho(psi[0]) ]

C = -args.alpha/args.epsilon**2
psi0 = [ args.Psi, -sqrt(2*rhoC*int_rho(args.Psi)), 0, 0 ]
sol = solve_ivp(f, [0, np.infty], psi0, dense_output=True)

X = np.linspace(0, sol.t[-1], args.N)
Psi = sol.sol(X)[0]
fig, axs = plt.subplots(ncols=2, figsize=args.figsize[0]*np.array((2.2,1)))
axs[0].plot(sol.t, sol.y[0], '*')
axs[0].plot(sol.t, sol.y[3], '*')
axs[0].plot(X, Psi, color='C0')
axs[0].set_title('psi')
axs[0].axhline(y=0, **Style.thin)
axs[1].plot(sol.t, rho(sol.y[0]), '*')
axs[1].plot(X, rho(Psi), color='C0')
axs[1].set_title('rho')

if args.pdf:
    plt.savefig('profile.pdf', bbox_inches='tight')
else:
    plt.show()

