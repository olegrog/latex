#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3.5, 3.5]}
py.rcParams.update(params)
import numpy as np
import sys, os, re, math, vtk, argparse
from scipy.interpolate import griddata
from vtk.util.numpy_support import vtk_to_numpy

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
parser.add_argument('latex', help='a LaTeX formula of the field')
parser.add_argument('--factor', default='1', metavar='expr', help='a multiplier for the field')
parser.add_argument('--kn', type=float, default=1, metavar='value', help='a Knudsen number for transform (x,y)-coords')
parser.add_argument('--lpow', type=float, default=1, metavar='value', help='a power for the levels arrangement')
parser.add_argument('--lmax', type=float, default=0.1, metavar='value', help='a shift for the maximum level')
parser.add_argument('--lmin', type=float, default=0, metavar='value', help='a shift for the minimum level')
parser.add_argument('--grid', type=int, default=200, metavar='value', help='a number of points along each axis')
args = parser.parse_args()
uniform_grid = True
grayscale = False

### Read data from VTK-file
reader = {
    '.vtu': vtk.vtkXMLUnstructuredGridReader,
    '.vtk': vtk.vtkUnstructuredGridReader
}[os.path.splitext(args.vtkfile)[1]]()
reader.SetFileName(args.vtkfile)
reader.Update()
x, y = get_points(reader)
F = get_point_data(reader, args.field) * eval(args.factor)

### Reconstruct the rectilinear grid from the unstructured one
xmin, xmax = min(x), max(x)
ymin, ymax = min(y), max(y)
if uniform_grid:
    N = args.grid
    xi = np.linspace(xmin, xmax, N)
    yi = np.linspace(ymin, ymax, N)
else:
    xi, yi = np.unique(x), np.unique(y)
    xi, yi = append_midpoints(xi), append_midpoints(yi)

### Draw scalar field
if len(F.shape) == 1:
    Fi = griddata((x, y), F, (xi[None,:], yi[:,None]), method='cubic')
    field_name = re.sub('Mean$', '', re.findall(r'[a-zA-Z]+', args.field)[0])
    lmin, lmax, laux = math.floor(10*np.min(Fi))/10, (math.floor(10*np.max(Fi)) + 1)/10, []
    lsteps = int(round(10*(lmax-lmin)))
    labelpos = False
    if field_name == 'T':
        laux = [ 1.05 ]
        R, npts = 0.45*xmax, lsteps + len(laux) - 1
        phi_min, phi_max = args.lmin, np.pi + args.lmax
        lx = xmax/2 + R*np.cos(np.linspace(phi_min, phi_max, npts))
        ly = R*np.sin(np.linspace(phi_min, phi_max, npts))
        labelpos = [(lx[i], ly[i]) for i in xrange(0, npts)]
    linewidth, fontsize, detailed = 1, 8, 4

### Draw vector field
else:
    Ui = griddata((x, y), F[:,0], (xi[None,:], yi[:,None]), method='cubic')
    Vi = griddata((x, y), F[:,1], (xi[None,:], yi[:,None]), method='cubic')
    Fi = np.sqrt(Ui*Ui + Vi*Vi)

    factor = 10**math.floor(np.log10(np.max(Fi)))
    lmin, laux = 0, []
    lmax = (math.floor(np.max(Fi)/factor*10) + 1.)*factor/10
    lsteps = int(round(10*lmax/factor))
    if lsteps >= 15 and lsteps < 30:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/2.))
        lmax *= 2.*lsteps/old_lsteps
    if lsteps >= 30 and lsteps < 75:
        old_lsteps = lsteps
        lsteps = int(round(lsteps/5.)) + 1
        lmax *= 5.*lsteps/old_lsteps
    if lsteps >= 75:
        old_lsteps = lsteps
        lsteps = int(math.floor(lsteps/10.)) + 1
        lmax *= 10.*lsteps/old_lsteps

    ### Calculate label positions along a line from min to max
    X, Y = np.meshgrid(xi, yi)
    area_min = np.where((np.abs(X-xmax/2) < xmax/4) & (Y < ymax/2) & (Y > ymax/8))
    area_max = np.where((np.abs(X-xmax/2) < xmax/4) & (Y < ymax/2))
    amin, amax = np.argmin(Fi[area_min]), np.argmax(Fi[area_max])
    Xmin, Ymin = X[area_min].flat[amin], Y[area_min].flat[amin]
    Xmax, Ymax = X[area_max].flat[amax], Y[area_max].flat[amax]
    lsteps_ = lsteps if np.sqrt((Xmax-Xmin)**2 + (Ymax-Ymin)**2) > np.sqrt((ymax-xmin)**2 + (ymax-ymin)**2)/2 else int(lsteps/2)+1
    lx = np.linspace(Xmin, Xmax, 1 + lsteps_)
    ly = np.linspace(Ymin**(1./args.lpow), Ymax**(1./args.lpow), 1 + lsteps_)**args.lpow
    #py.plot(lx, ly, 'D')
    labelpos = [(lx[i], ly[i]) for i in xrange(1, lsteps_)]
    linewidth, fontsize, detailed = .5, 6, 2


print '%s: max(magF) = %g, lmin = %g, lmax = %g, lsteps = %g' % \
    (args.pdffile, np.max(Fi), lmin, lmax, lsteps)

### Plot color detailed contours
levels = np.linspace(lmin, lmax, 1 + lsteps*detailed)
if grayscale:
    py.contour(xi, yi, Fi, levels=levels, colors='lightgray', linewidths=0.05)
else:
    cmap = py.cm.get_cmap('coolwarm')
    py.contourf(xi, yi, Fi, levels=levels, cmap=cmap)

### Plot black contours with labels
levels = np.linspace(lmin, lmax, 1 + lsteps)
levels = np.sort(np.append(levels, laux))
CS = py.contour(xi, yi, Fi, levels=levels, colors='k', linewidths=linewidth)
clabels = py.clabel(CS, levels,
    inline=True,
    use_clabeltext=True,
    fmt='%g',
    manual=labelpos,
    fontsize=fontsize)
[ txt.set_backgroundcolor('white') for txt in clabels ]
[ txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels ]

### Add a streamplot for the vector field
# NB: py.streamplot relies on evenly grid
if len(F.shape) > 1:
    py.streamplot(xi, yi, Ui, Vi,
        density=0.89,
        minlength=.2,
        arrowstyle='->',
        color='k',
        linewidth=1)

py.text(-.05, .4, r'$\displaystyle %s$' % args.latex)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$y$', labelpad=-5, rotation=0)
py.xlim(xmin, xmax)
py.ylim(ymin, ymax)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([xmin, xmax])
ax.set_xticks([ymin, ymax])
py.tick_params(axis='both', direction='out')
py.savefig(args.pdffile, bbox_inches='tight', transparent=True)
