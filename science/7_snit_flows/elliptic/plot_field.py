#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'font.size': 10,
          'axes.labelsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [7.5,4]}
py.rcParams.update(params)
import numpy as np
import sys, os, re, math, vtk, argparse
from datetime import datetime
from matplotlib.tri import Triangulation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse, Arc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter

parser = argparse.ArgumentParser(description='2D contour plots of scalar and vector from VTK data')
parser.add_argument('vtkfile', help='an input VTK file')
parser.add_argument('pdffile', help='an output PDF file')
parser.add_argument('field', help='a name of the VTK field')
parser.add_argument('--latex', default='\\frac{v_i}{k}', metavar='latex', help='a LaTeX formula of the field')
parser.add_argument('--factor', default='1', metavar='expr', help='a multiplier for the field')
parser.add_argument('--kn', type=float, default=1, metavar='value', help='a Knudsen number for transform (x,y)-coords')
parser.add_argument('--lmax', type=float, default=1e6, metavar='value', help='maximum contour level')
parser.add_argument('--lpow', type=float, default=1, metavar='value', help='a power for the levels arrangement')
args = parser.parse_args()
grayscale = True

a1 = 1.5
a0 = 0.3
b0 = 0.7

def get_point_data(reader, field):
    vtk_array = reader.GetOutput().GetPointData().GetArray(field)
    return vtk_to_numpy(vtk_array)

def get_triangles(reader):
    vtk_array = reader.GetOutput().GetPolys().GetData()
    N = reader.GetOutput().GetNumberOfPolys()
    return np.delete(vtk_to_numpy(vtk_array).reshape((N, 4)), 0, 1)

def get_points(reader):
    vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array) * args.kn
    return numpy_array[:,0], numpy_array[:,1], numpy_array[:,2]

def get_levels(maxU):
    factor = 10**math.floor(np.log10(maxU))
    lmin, laux = 0, []
    lmax = min((math.floor(maxU/factor*10) + 1.)*factor/10, args.lmax)
    lsteps = int(round(10*lmax/factor))
    if lsteps >= 15 and lsteps < 30:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/2.))
        lmax *= 2.*lsteps/old_lsteps
    if lsteps >= 30 and lsteps < 60:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/5.)) + 1
        lmax *= 5.*lsteps/old_lsteps
    if lsteps >= 60:
        old_lsteps = lsteps
        lsteps = int(math.floor(lsteps/10.)) + 1
        lmax *= 10.*lsteps/old_lsteps
    return lmin, lmax, lsteps

def draw_artists(ax, beta, lw):
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

### Read data from VTK-file
reader = {
    '.vtu': vtk.vtkXMLUnstructuredGridReader,
    '.vtk': vtk.vtkUnstructuredGridReader
}[os.path.splitext(args.vtkfile)[1]]()
reader.SetFileName(args.vtkfile)
reader.Update()

### Create a plane cut with triangulation
plane = vtk.vtkPlane()
plane.SetOrigin(reader.GetOutput().GetCenter())
plane.SetNormal(0, 0, 1)
planeCut = vtk.vtkCutter()
planeCut.SetInputConnection(reader.GetOutputPort())
planeCut.SetCutFunction(plane)
planeCut.GenerateTrianglesOn()      # return triangles instead of quads
planeCut.Update()

### Get data
X, Y, Z = get_points(planeCut)
#F = gaussian_filter(get_point_data(planeCut, args.field), (1,0)) * eval(args.factor)
F = get_point_data(planeCut, args.field) * eval(args.factor)
tris = get_triangles(planeCut)
U, V = F[:,0], F[:,1]
magU = np.sqrt(U*U+V*V)
lmin, lmax, lsteps = get_levels(np.max(magU))
triang = Triangulation(X, Y, tris)
print '%s: max(magF) = %g, lmin = %g, lmax = %g, lsteps = %g' % \
    (args.pdffile, np.max(magU), lmin, lmax, lsteps)

### Calculate label positions along a line from min to max
xmin, xmax = min(X), max(X)
ymin, ymax = min(Y), max(Y)
area_min = np.where((X > xmax/2) & (X < .8*xmax) & (Y < ymax/4))
area_max = np.where((X+Y > 0.9))
amin, amax = np.argmin(magU[area_min]), np.argmax(magU[area_max])
Xmin, Ymin = X[area_min].flat[amin], Y[area_min].flat[amin]
Xmax, Ymax = X[area_max].flat[amax], Y[area_max].flat[amax]
lsteps_ = int(max(magU[area_max])/lmax*lsteps) + 1
lx = np.linspace(Xmin, Xmax, 1 + lsteps_)
ly = np.linspace(Ymin**(1./args.lpow), Ymax**(1./args.lpow), 1 + lsteps_)**args.lpow
#py.plot(lx, ly, 'D')
labelpos = [(lx[i], ly[i]) for i in xrange(1, lsteps_)]
linewidth, fontsize = .5, 6

### Plot detailed contours
levels = np.linspace(lmin, lmax, 1 + lsteps*2)
if grayscale:
    py.tricontour(triang, magU, levels=levels, colors='lightgray', linewidths=0.05)
else:
    cmap = py.cm.get_cmap('coolwarm')
    py.tricontourf(triang, magU, cmap=cmap, levels=levels)

### Plot black contours with labels
levels = np.linspace(lmin, lmax, 1 + lsteps)
CS = py.tricontour(triang, magU, levels=levels, colors='k', linewidths=.3)
clabels = py.clabel(CS, levels,
    inline=True,
    use_clabeltext=True,
    fmt='%g',
    manual=labelpos,
    fontsize=6)
[ txt.set_backgroundcolor('white') for txt in clabels ]
[ txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.5)) for txt in clabels ]
### Draw a streamplot
# NB: py.streamplot relies on evenly grid => create it!
N, method = 100, 'cubic'
xi = np.linspace(X.min(), X.max(), (X.max()-X.min())*N)
yi = np.linspace(Y.min(), Y.max(), (Y.max()-Y.min())*N)
print "before griddata", datetime.now().time()
U = griddata((X, Y), U, (xi[None,:], yi[:,None]), method=method)
V = griddata((X, Y), V, (xi[None,:], yi[:,None]), method=method)
print "after griddata", datetime.now().time()
magU = np.sqrt(U*U+V*V)
lw = 1 #5*magU/magU.max()
kwargs = {}
if grayscale:
    kwargs['color'] = 'k'
else:
    kwargs = {
        'color': magU/magU.max(),
        'cmap': py.cm.get_cmap('autumn')
    }
py.streamplot(xi, yi, U, V,
    density=0.9,
    minlength=.2,
    arrowstyle='->',
    linewidth=lw,
    **kwargs)

### Draw additional objects
ax = py.axes()
draw_artists(ax, np.pi/2, 0.75)

py.setp(ax,
    xticks=[0, a0, a1], xticklabels=[r'$0$', r'$a_0$', r'$a_1$'],
    yticks=[0, b0, 1], yticklabels=[r'$0$', r'$b_0$', r'$1$'])

ax.text(1.5, .9, r'$\displaystyle %s$' % args.latex)
ax.text(.15, .4, r'$T_1$', zorder=50)
ax.text(1., .8, r'$T_2$')
ax.set_xlim(0, a1+.05)
ax.set_ylim(0, 1+.05)
ax.text(a1/2, -.09, r'$x$')
ax.text(-.09, .85, r'$y$')

py.tick_params(axis='both', direction='out',)
py.savefig(args.pdffile, bbox_inches='tight', transparent=True)
