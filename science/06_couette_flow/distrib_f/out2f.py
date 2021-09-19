#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 9,
          'font.size': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'figure.figsize': [3.0,2.2]}
py.rcParams.update(params)
import numpy as np
import sys, math, array, struct, itertools, json
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline

sys.path.append('/Users/olegrog/kesolver/tools')
from kepy.ximesh import read_ximesh
from out2 import readNodesCells, center

keifile = sys.argv[1]
binfile = sys.argv[2]
outfile = sys.argv[3]
#celli = int(sys.argv[4])
position = sys.argv[4]
smooth = float(sys.argv[5]) if len(sys.argv) > 5 else 0
axis = sys.argv[6] if len(sys.argv) > 6 else 'z'
node = int(sys.argv[7]) if len(sys.argv) > 7 else 0
cut2 = int(sys.argv[8]) if len(sys.argv) > 8 else 3.0

def slice2d(f, ax, rad):
    rad += node
    return {
        'x': f[rad,:,:].T,
        'y': f[:,rad,:].T,
        'z': f[:,:,rad]
    }[ax]

def slice1d(f, g, ax, rad):
    idx = rad + node
    return {
        'x': (f[idx,rad,:], g[idx,:,rad]),
        'y': (f[rad,idx,:], g[:,idx,rad]),
        'z': (f[:,rad,idx], g[rad,:,idx])
    }[ax]

def axis_pair(ax):
    return {
        'x': ('z','y'),
        'y': ('z','x'),
        'z': ('x','y')
    }[ax]

ax2int = lambda ax: ord(ax) - ord('x')

ax = py.axes()
with open(keifile, 'rb') as fd:
    data = json.load(fd)
nodes, cells = readNodesCells(keifile)
N = len(cells)
symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(data)
xyz /= np.sqrt(2)
vol /= 2*np.sqrt(2)
size = np.sum(circl)
f = np.zeros_like(r)
f.fill(np.NaN)
f[circl] = 0
our = axis_pair(axis)

ycoords = [center([nodes[n] for n in cells[i].nodes])[1] for i in range(N)]
ycoord = {
    'center': np.min(np.abs(ycoords)),
    'boundary': np.max(ycoords)
}[position]
celli = np.argwhere(np.abs(ycoords) == ycoord).flatten()
spline_type = 'cubic' if smooth > 0 else 'linear'
k_spline = {
    'cubic': 3,
    'linear': 1
}[spline_type]

with open(binfile, 'rb') as fd:
    for _ in range(N):
        data = fd.read(4+4)
        i, _ = struct.unpack('=ii', data)
        a = array.array('d')
        a.fromfile(fd, size)
        if i in celli:
            if center([nodes[n] for n in cells[i].nodes])[1] > 0:
                print('up', i, np.sum(a))
                f[circl] += np.array(a)
            else:
                print('down', i, np.sum(a))
                f[::-1,::-1,:][circl] += np.array(a)
            
f /= len(celli)
idx, d = [rad,rad,rad], ax2int(axis)
idx[d] += node
print("Total mass = %f, y_coord = %f, filename = %s, zeta_%s = %f, smooth = %.2e" % \
    (np.nansum(f), ycoord, outfile, axis, xyz[d][tuple(idx)], smooth), file=sys.stderr)
f /= vol

gap = 0 #0.3
cut = cut2 - gap
mesh_size = lambda ax: 151 if ax == 'y' else 91
ratio = lambda ax: 1.8 if ax == 'y' else 1
base_x, base_y = np.linspace(-1, 1, mesh_size(our[0])), np.linspace(-1, 1, mesh_size(our[1]))
# [ py.axhline(y, c='k', lw=0.05) for y in np.sign(base_y)*np.abs(base_y)**ratio(our[1])*cut ]
# [ py.axvline(x, c='k', lw=0.05) for x in np.sign(base_x)*np.abs(base_x)**ratio(our[0])*cut ]
grid_x, grid_y = np.meshgrid(
    np.sign(base_x)*np.abs(base_x)**ratio(our[0])*cut,
    np.sign(base_y)*np.abs(base_y)**ratio(our[1])*cut,
    sparse=False, indexing='ij')
X, Y = slice1d(xyz[ax2int(our[0])], xyz[ax2int(our[1])], axis, rad)
fnew = np.nan_to_num(f)
mX, mY = np.abs(X) < cut, np.abs(Y) < cut
X, Y = X[mX], Y[mY]
F = slice2d(fnew, axis, rad)[np.outer(mX,mY)].reshape(len(X), len(Y))
spline = RectBivariateSpline(X, Y, F, kx=k_spline, ky=k_spline, s=smooth)
grid_z = np.zeros_like(grid_x)
grid_z.fill(np.NaN)
mask = grid_x**2 + grid_y**2 < 10*cut**2
grid_z[mask] = np.maximum(spline(grid_x[mask], grid_y[mask], grid=False), 0)
#zmax = np.round(np.nanmax(grid_z)*1.01, 1-int(np.log(np.nanmax(grid_z))))
rnd = 500
zmax = 1./rnd*(int(rnd*np.nanmax(grid_z))+1)
print("zmin = %.3f, zmax = %.3f, round = %.3f" % \
    (np.nanmin(grid_z), np.nanmax(grid_z), zmax), file=sys.stderr)
ax.set_xlim(-cut2, cut2)
ax.set_ylim(-cut2, cut2)
levels = np.linspace(0, zmax, 21)
cset = ax.contourf(grid_x, grid_y, grid_z, zdir='z', cmap=py.cm.get_cmap('coolwarm'), levels=levels)
cbar = py.colorbar(cset, ticks=levels[::4], drawedges=False)
cbar.solids.set_edgecolor("face")
ax.set_xlabel(r'$\xi_%s$' % our[0])
ax.set_ylabel(r'$\xi_%s$' % our[1], rotation=0)
ax.set_aspect('equal')
dashes = [1,1]
py.axhline(0, c='k', lw=0.25, dashes=dashes)
py.axvline(0, c='k', lw=0.25, dashes=dashes)

py.savefig(outfile, bbox_inches='tight', pad_inches=0.01, transparent=True)
#py.show()

