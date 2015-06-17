#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [6,5]}
py.rcParams.update(params)
import numpy as np
import sys, os, re, math, vtk
from scipy.interpolate import griddata
from vtk.util.numpy_support import vtk_to_numpy

infile = sys.argv[1]
outfile = sys.argv[2]
field = sys.argv[3]
tex_name = sys.argv[4]
kn = float(sys.argv[5]) if len(sys.argv) > 5 else -1
kpow = float(sys.argv[6]) if len(sys.argv) > 6 else 0
lpow = float(sys.argv[7]) if len(sys.argv) > 7 else 1

uniform_grid = True

def get_point_data(reader, field):
    vtk_array = reader.GetOutput().GetPointData().GetArray(field)
    return vtk_to_numpy(vtk_array)

def get_points(reader):
    vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array)
    return numpy_array[:,0], numpy_array[:,1]

def append_midpoints(x):
    return np.sort(np.append(x, .5*(x[1:] + x[:-1])))

### Read data from VTK-file
reader = {
    '.vtu': vtk.vtkXMLUnstructuredGridReader,
    '.vtk': vtk.vtkUnstructuredGridReader
}[os.path.splitext(infile)[1]]()
reader.SetFileName(infile)
reader.Update()
x, y = get_points(reader)
F = get_point_data(reader, field)

### Some corrections to data
if kn > 0:
    k = kn*np.sqrt(np.pi)/2
    F *= k**kpow
if outfile.startswith('data'):
    x, y = kn*x, kn*y 
    y = 0.5 - y

### Reconstruct the rectilinear grid from the unstructured one
xmin, xmax = min(x), max(x)
ymin, ymax = min(y), max(y)
if uniform_grid:
    N = 200
    xi = np.linspace(xmin, xmax, N)
    yi = np.linspace(ymin, ymax, N)
else:
    xi, yi = np.unique(x), np.unique(y)
    xi, yi = append_midpoints(xi), append_midpoints(yi)

### Draw scalar field
if len(F.shape) == 1:
    Fi = griddata((x, y), F, (xi[None,:], yi[:,None]), method='cubic')
    field_name = re.findall(r'[a-zA-Z]+', field)[0]
    lmin, lmax, laux = {
        'T': (.5, 1.5, [ 1.05 ]),
        'rho': (.5, 2., [])
    }[field_name]
    lsteps = int(round(10*(lmax-lmin)))
    labelpos = False
    if field_name == 'T':
        R, npts = 0.45*xmax, lsteps + len(laux) - 1
        phi_min, phi_max, r = -0.05, np.pi + 0.15, 0
        if np.max(Fi) < 1.4:
            r = 1
        lx = xmax/2 + R*np.cos(np.linspace(phi_min, phi_max, npts))
        ly = R*np.sin(np.linspace(phi_min, phi_max, npts))
        labelpos = [(lx[i], ly[i]) for i in xrange(0+r, npts-r)]
    linewidth, fontsize = 1, 8

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

    ### Calculate label positions along a line from min to max
    X, Y = np.meshgrid(xi, yi)
    area = np.where((np.abs(X-xmax/2) < xmax/4) & (Y < ymax/2))
    amin, amax = np.argmin(Fi[area]), np.argmax(Fi[area])
    X, Y = X[area].flat, Y[area].flat
    lx = np.linspace(X[amin], X[amax], 1 + lsteps)
    ly = np.linspace(Y[amin]**(1./lpow), Y[amax]**(1./lpow), 1 + lsteps)**lpow
    labelpos = [(lx[i], ly[i]) for i in xrange(1, lsteps)]
    linewidth, fontsize = .5, 6


print '%s: max(magF) = %g, lmin = %g, lmax = %g, lsteps = %g' % \
    (outfile, np.max(Fi), lmin, lmax, lsteps)

### Plot color detailed contours
cmap = py.cm.get_cmap('coolwarm')
levels = np.linspace(lmin, lmax, 1 + lsteps*8)
py.contourf(xi, yi, Fi, levels=levels, cmap=cmap)

### Plot black contours with labels
levels = np.linspace(lmin, lmax, 1 + lsteps)
levels = np.append(levels, laux)
CS = py.contour(xi, yi, Fi, levels=levels, colors='k', linewidths=linewidth)
py.clabel(CS, levels,
    inline=True,
    use_clabeltext=True,
    fmt='%g',
    manual=labelpos,
    fontsize=fontsize)

### Add a streamplot for the vector field
# NB: py.streamplot relies on evenly grid
if len(F.shape) > 1:
    lw = 1 #5*Fi/Fi.max()
    py.streamplot(xi, yi, Ui, Vi,
        color=Fi/Fi.max(),
        density=0.89,
        minlength=.2,
        arrowstyle='->',
        linewidth=lw,
        cmap=py.cm.get_cmap('autumn'))

py.text(-.05, .4, r'$\displaystyle %s$' % tex_name)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$y$', labelpad=-5, rotation=0)
py.xlim(xmin, xmax)
py.ylim(ymin, ymax)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([xmin, xmax])
ax.set_xticks([ymin, ymax])
py.tick_params(axis='both', direction='out')
py.savefig(outfile, bbox_inches='tight', transparent=True)
