#!/usr/bin/env python

import matplotlib, scipy.ndimage, argparse
import numpy as np
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'backend': 'pdf',
    'font.size': 10,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.linewidth': 1,
    'text.usetex': True,
    'figure.figsize': [4,4]
})

parser = argparse.ArgumentParser(description='2D interpolated contour plots')
parser.add_argument('txtfile', help='an input txt file')
parser.add_argument('pdffile', help='an output PDF file')
parser.add_argument('Nx', type=int, help='a number of points along the x axis')
parser.add_argument('Ny', type=int, help='a number of points along the y axis')
parser.add_argument('-x', default='x', help='formula for the x axis')
parser.add_argument('-y', default='y', help='formula for the y axis')
parser.add_argument('--zoom', type=int, default=3, help='zoom factor')
parser.add_argument('--column', type=int, default=3, help='column of the values')
parser.add_argument('--ground', type=float, default=0, help='zero level for the function')
parser.add_argument('--xlabel', default='$x$', help='label for the x axis')
parser.add_argument('--ylabel', default='$y$', help='label for the y axis')
parser.add_argument('--logx', action='store_true', help='log scaling on the x axis')
parser.add_argument('--logy', action='store_true', help='log scaling on the y axis')
args = parser.parse_args()

refine = lambda A: scipy.ndimage.zoom(np.reshape(A, (args.Ny, args.Nx)), args.zoom)

data = np.loadtxt(args.txtfile)
data = data[data[:,1].argsort()]
data = data[data[:,0].argsort(kind='mergesort')].T
x, y, f = data[0], data[1], data[args.column-1] - args.ground
x_ = eval(args.x); y_ = eval(args.y)
x = x_; y = y_
plt.scatter(x, y, marker='o', s=5, c='k', zorder=10)

x, y, f = refine(x), refine(y), refine(f)
plt.contourf(x, y, f, levels=50, cmap = plt.cm.get_cmap('coolwarm'))
CS = plt.contour(x, y, f, levels=10, colors='k')
plt.clabel(CS, inline=True, fmt='%g', fontsize=10)

#ax = plt.gca()
#ax.set_aspect(1./ax.get_data_ratio())
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(y), np.max(y))

plt.xlabel(args.xlabel)
plt.ylabel(args.ylabel)

if args.logx:
    plt.semilogx()
if args.logy:
    plt.semilogy()
plt.tick_params(axis='both', direction='out',)
plt.savefig(args.pdffile, bbox_inches='tight', transparent=True)
