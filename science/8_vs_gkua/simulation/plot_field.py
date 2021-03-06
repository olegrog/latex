#!/usr/bin/env python

import sys, os, re, math, vtk, argparse, subprocess
import numpy as np
import matplotlib as mpl
from matplotlib import tri, cm, pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata

params = {'backend': 'pdf',
          'font.size': 10,
          'axes.labelsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [7.5,4]}
mpl.rcParams.update(params)

vec = lambda s: np.fromstring(s, sep=',')
parser = argparse.ArgumentParser(description='2D contour plots of scalar and vector from VTK data')
parser.add_argument('vtkfile', help='an input VTK file')
parser.add_argument('pdffile', help='an output PDF file')
parser.add_argument('field', help='a name of the VTK field')
parser.add_argument('--latex', default='\\mathrm{Ma}', metavar='latex', help='a LaTeX formula of the field')
parser.add_argument('--factor', default='1', metavar='expression', help='a multiplier for the field')
parser.add_argument('--refine', type=int, default='2', metavar='value', help='each triangle is divided into 4**n child triangles')
parser.add_argument('--lmax', type=float, default=1e6, metavar='value', help='maximum contour level')
parser.add_argument('--lsteps-max', type=int, default=10, metavar='value', help='maximum number of contour levels')
parser.add_argument('--lpow', type=float, default=1, metavar='value', help='a power for the levels arrangement')
parser.add_argument('--lrad', type=float, default=0.1, metavar='value', help='radius of searching extrema')
parser.add_argument('--lcmax', type=vec, default='0,0', metavar='vector', help='coordinates of maximum')
parser.add_argument('--lcmin', type=vec, default='0,0', metavar='vector', help='coordinates of minimum')
parser.add_argument('--show', action='store_true', help='if label points are shown')
parser.add_argument('--draw-mesh', action='store_true', help='if mesh is drawn')
parser.add_argument('--white-clab', action='store_true', help='if clabels is placed in the white box')
parser.add_argument('--export-dat', action='store_true', help='export to Tecplot dat-file instead of plotting')
args = parser.parse_args()
grayscale = False

class geom:
    r = eval(subprocess.check_output('grep "^r =" $(dirname %s)/../*.geo | awk \'{ print $3 }\' | sed "s/;//"' % args.vtkfile, shell=True), { 'Kn': 1 })
    R = eval(subprocess.check_output('grep "^R =" $(dirname %s)/../*.geo | awk \'{ print $3 }\' | sed "s/;//"' % args.vtkfile, shell=True), { 'r': r })

def get_point_data(reader, field):
    vtk_array = reader.GetOutput().GetPointData().GetArray(field)
    return vtk_to_numpy(vtk_array)

def get_triangles(reader):
    vtk_array = reader.GetOutput().GetPolys().GetData()
    N = reader.GetOutput().GetNumberOfPolys()
    return np.delete(vtk_to_numpy(vtk_array).reshape((N, 4)), 0, 1)

def get_points(reader):
    vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array)
    return numpy_array[:,0], numpy_array[:,1], numpy_array[:,2]

def get_levels(maxU):
    factor = 10**math.floor(np.log10(maxU))
    lmin, laux = 0, []
    lmax = min((int(maxU/factor*10) + 1.)*factor/10, args.lmax)
    lsteps = int(round(10*lmax/factor))
    if lsteps >= args.lsteps_max and lsteps < 2*args.lsteps_max:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/2.))
        lmax *= 2.*lsteps/old_lsteps
    if lsteps >= 2*args.lsteps_max and lsteps < 5*args.lsteps_max:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/5.))
        lmax *= 5.*lsteps/old_lsteps
    if lsteps >= 5*args.lsteps_max and lsteps < 10*args.lsteps_max:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/10.))
        lmax *= 10.*lsteps/old_lsteps
    if lsteps >= 10*args.lsteps_max and lsteps < 20*args.lsteps_max:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/20.))
        lmax *= 20.*lsteps/old_lsteps
    return lmin, lmax, lsteps

def draw_artists(ax, beta, lw):
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    artists = [
        Circle(xy=[0,0], radius=geom.R, color='k', fill=False, lw=1.5*lw),
        Circle(xy=[0,0], radius=geom.r, color='w', fill=True, zorder=20),
        Circle(xy=[0,0], radius=geom.r, color='k', fill=False,lw=1.5*lw,  zorder=20),
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
plane.SetOrigin(0, 0, 1e-13) 
plane.SetNormal(0, 0, 1)
planeCut = vtk.vtkCutter()
planeCut.SetInputConnection(reader.GetOutputPort())
planeCut.SetCutFunction(plane)
planeCut.GenerateTrianglesOn()      # return triangles instead of quads
planeCut.Update()

### Get data
X, Y, Z = get_points(planeCut)
F = get_point_data(planeCut, args.field) * eval(args.factor)
tris = get_triangles(planeCut)
U, V = F[:,0], F[:,1]
magU = np.sqrt(U*U+V*V)
lmin, lmax, lsteps = get_levels(np.max(magU))
triang = tri.Triangulation(X, Y, tris)
#plt.triplot(triang, lw=0.1, color='red')

print '%s: max(magU) = %g, lmin = %g, lmax = %g, lsteps = %g' % \
    (args.pdffile, max(magU), lmin, lmax, lsteps)

### Refine triangulation and data
refiner = tri.UniformTriRefiner(triang)
# default kind='min_E' does not work for bad quality triangles
tri_interp = tri.CubicTriInterpolator(triang, magU, kind='geom')
refi_triang, refi_magU = refiner.refine_field(magU, subdiv=args.refine, triinterpolator=tri_interp)

### Calculate label positions along a line from min to max
xmin, xmax = min(X), max(X)
ymin, ymax = min(Y), max(Y)
area_min = np.where((X-args.lcmin[0])**2 + (Y-args.lcmin[1])**2 < args.lrad**2)
area_max = np.where((X-args.lcmax[0])**2 + (Y-args.lcmax[1])**2 < args.lrad**2)
amin, amax = np.argmin(refi_magU[area_min]), np.argmax(refi_magU[area_max])
Xmin, Ymin = X[area_min].flat[amin], Y[area_min].flat[amin]
Xmax, Ymax = X[area_max].flat[amax], Y[area_max].flat[amax]
lsteps_ = int(max(refi_magU[area_max])/lmax*lsteps)
lx = Xmin + np.sign(Xmax-Xmin)*np.linspace(0, abs(Xmax-Xmin)**(1./args.lpow), 1 + lsteps_)**args.lpow
ly = Ymin + np.sign(Ymax-Ymin)*np.linspace(0, abs(Ymax-Ymin)**(1./args.lpow), 1 + lsteps_)**args.lpow
if args.show:
    plt.plot(lx, ly, '*', ms=3, zorder=100)
lx, ly = (lx[1:] + lx[:-1])/2, (ly[1:] + ly[:-1])/2
if args.show:
    plt.plot(lx, ly, 'D', ms=3, zorder=100)
labelpos = tuple(map(tuple, np.vstack((lx,ly)).T))
linewidth, fontsize = .5, 6

### Plot detailed contours
levels = np.linspace(lmin, lmax, 1 + lsteps*2)
if grayscale:
    plt.tricontour(refi_triang, refi_magU, levels=levels, colors='lightgray', linewidths=0.05)
else:
    cmap = cm.get_cmap('coolwarm')
    plt.tricontourf(refi_triang, refi_magU, cmap=cmap, levels=levels)

### Plot black contours with labels
levels = np.linspace(lmin, lmax, 1 + lsteps)
CS = plt.tricontour(refi_triang, refi_magU, levels=levels, colors='k', linewidths=.3)
if args.draw_mesh:
    plt.triplot(refi_triang, lw=0.05, color='gray')
    plt.triplot(triang, lw=0.1, color='k')

clabels = plt.clabel(CS, levels,
    inline=True,
    use_clabeltext=True,
    fmt='%g',
    manual=labelpos,
    fontsize=5, #6
    zorder=50)
if args.white_clab:
    [ txt.set_backgroundcolor('white') for txt in clabels ]
    [ txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.5)) for txt in clabels ]

### Draw a streamplot
# NB: plt.streamplot relies on evenly grid => create it!
N, method = 50, 'cubic'
xi = np.linspace(xmin, xmax, (xmax-xmin)*N)
yi = np.linspace(ymin, ymax, (ymax-ymin)*N)
U = griddata((X, Y), U, (xi[None,:], yi[:,None]), method=method)
V = griddata((X, Y), V, (xi[None,:], yi[:,None]), method=method)

magU = np.sqrt(U*U+V*V)
lw = 1 #5*magU/magU.max()
kwargs = {}
if grayscale:
    kwargs['color'] = 'k'
else:
    kwargs = {
        'color': magU/magU.max(),
        'cmap': cm.get_cmap('autumn')
    }
plt.streamplot(xi, yi, U, V,
    density=0.5,
    minlength=.2,
    arrowstyle='->',
    linewidth=lw*.5)
#    ,
#    **kwargs)

### Draw additional objects
ax = plt.axes()
draw_artists(ax, np.pi/2, 0.75)

plt.setp(ax,
    xticks=[-geom.R, -geom.r, geom.r, geom.R], xticklabels=[r'$-R$', r'$-r$', r'$r$', r'$R$'])

ax.set_xlim(-geom.R-.5, geom.R+.05)
ax.set_ylim(0, geom.R+.05)
ax.text(geom.R-.5, geom.R-.5, r'$\displaystyle %s$' % args.latex)

plt.tick_params(axis='both', direction='out',)
plt.savefig(args.pdffile, bbox_inches='tight', transparent=True)
