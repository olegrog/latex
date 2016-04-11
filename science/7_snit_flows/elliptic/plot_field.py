#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'font.size': 10,
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [7.5,4]}
py.rcParams.update(params)
import numpy as np
import sys, os, re, math, vtk, argparse
from matplotlib.tri import Triangulation
from matplotlib.mlab import griddata
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse, Arc
from vtk.util.numpy_support import vtk_to_numpy

a1 = 1.5
a0 = 0.3
b0 = 0.7

def get_point_data(reader, field):
    vtk_array = reader.GetOutput().GetPointData().GetArray(field)
    return vtk_to_numpy(vtk_array)

def get_points(reader):
    vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array) * args.kn
    return numpy_array[:,0], numpy_array[:,1]

def append_midpoints(x):
    return np.sort(np.append(x, .5*(x[1:] + x[:-1])))

parser = argparse.ArgumentParser(description='2D contour plots of scalar and vector from VTK data')
parser.add_argument('vtkfile', help='an input VTK file')
parser.add_argument('pdffile', help='an output PDF file')
parser.add_argument('field', help='a name of the VTK field')
parser.add_argument('--factor', default='1', metavar='expr', help='a multiplier for the field')
parser.add_argument('--kn', type=float, default=1, metavar='value', help='a Knudsen number for transform (x,y)-coords')
args = parser.parse_args()
grayscale = False

### Read data from VTK-file
reader = {
    '.vtu': vtk.vtkXMLUnstructuredGridReader,
    '.vtk': vtk.vtkUnstructuredGridReader
}[os.path.splitext(args.vtkfile)[1]]()
reader.SetFileName(args.vtkfile)
reader.Update()
X, Y = get_points(reader)
F = get_point_data(reader, args.field) * eval(args.factor)

def draw(ax, beta, lw):
    U, V = F[:,0], F[:,1]
    magU = np.sqrt(U*U+V*V)

    triang = Triangulation(X,Y)
    levels = np.linspace(0, (int(magU.max()*100) + 1.)/100, 11)
    Fi = np.sqrt(U*U+V*V)
    if grayscale:
        CF = ax.tricontourf(triang, magU, levels=levels, colors='lightgray', linewidths=0.05)
    else:
        cmap = py.cm.get_cmap('coolwarm')
        CF = ax.tricontourf(triang, magU, cmap=cmap, levels=levels)

    N = 100
    interp = 'linear'
    xi = np.linspace(X.min(), X.max(), (X.max()-X.min())*N)
    yi = np.linspace(Y.min(), Y.max(), (Y.max()-Y.min())*N)
    U = griddata(X,Y,U,xi,yi,interp=interp)
    V = griddata(X,Y,V,xi,yi,interp=interp) 
    ax.streamplot(xi, yi, U, V, color='k', density=2*np.sqrt(Y.max()-Y.min()), minlength=.2, linewidth=lw, arrowstyle='->', arrowsize=lw)
    
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    angle = -(np.pi/2-beta)/np.pi*180
    artists = [
        Ellipse(xy=[0,0], width=2*a1, height=2, angle=0, color='k', fill=False, lw=1.5*lw),
        Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=angle, color='w', fill=True, zorder=20),
        Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=angle, color='k', fill=False, lw=1.5*lw, zorder=30),
    ]

    for a in artists:
        ax.add_artist(a)
    
    return CF

ax = py.axes()
CF = draw(ax, np.pi/2, 0.75)
cax = make_axes_locatable(ax).append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax)

py.setp(ax,
    xticks=[0, a0, a1], xticklabels=[r'$0$', r'$a_0$', r'$a_1$'],
    yticks=[0, b0, 1], yticklabels=[r'$0$', r'$b_0$', r'$1$'])

ax.text(1.5, .9, r'$\displaystyle\frac{v_i}{k}$')
ax.text(.15, .4, r'$T_1$', zorder=50)
ax.text(1., .8, r'$T_2$')
ax.set_xlim(0, a1+.05)
ax.set_ylim(0, 1+.05)
ax.text(a1/2, -.09, r'$x$')
ax.text(-.09, .85, r'$y$')

py.tick_params(axis='both', direction='out',)
py.savefig(args.pdffile, bbox_inches='tight', transparent=True)
