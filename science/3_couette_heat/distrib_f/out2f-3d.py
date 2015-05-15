#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import numpy as np

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline

sys.path.append('/Users/olegrog/kesolver/tools')
from kepy.ximesh import read_ximesh
from out2 import readNodesCells, center

keifile = sys.argv[1]
binfile = sys.argv[2]
outfile = sys.argv[3]
cell = int(sys.argv[4])
smooth = float(sys.argv[5]) if len(sys.argv) > 5 else 0
axis = sys.argv[6] if len(sys.argv) > 6 else 'z'
node = int(sys.argv[7]) if len(sys.argv) > 7 else 0
cut2 = int(sys.argv[8]) if len(sys.argv) > 8 else 4.0

def sqr(x):
    return x*x

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


def init_plot_params(columnwidth, height):
    fig_width_pt = columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = height*inches_per_pt
    fig_size =  [fig_width,fig_height]
    params = {
                'backend': 'pdf',
                'axes.labelsize': 10,
                'font.size': 10,
                'legend.fontsize': 9,
                'xtick.labelsize': 7,
                'ytick.labelsize': 7,
                'text.usetex': False,
                'figure.figsize': fig_size,
                'font.family':'serif'
            }
    plt.rcParams.update(params)

width = 350
height = width * 0.75
init_plot_params(width, height)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=-0.15, right=1.1, bottom=0, top=1.2)
with open(keifile, 'rb') as fd:
    data = json.load(fd)

nodes, cells = readNodesCells(keifile)

symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(data)
xyz /= np.sqrt(2)

with open(binfile, 'rb') as fd:
    while fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)
        f = np.zeros_like(r)
        f.fill(np.NaN) #-1e-10)
        a = array.array('d')
        size = np.sum(circl)
        a.fromfile(fd, size)

        if i == int(cell):
            f[circl] = np.array(a)
            idx, d = [rad,rad,rad], ax2int(axis)
            idx[d] += node
            print >> sys.stderr, "Total mass = %f, y_coord = %f, filename = %s, zeta_%s = %f" % \
                (np.nansum(f), center([nodes[n] for n in cells[i].nodes])[1], outfile, axis, xyz[d][tuple(idx)]/np.sqrt(2))
            f /= vol
            our = axis_pair(axis)
            xp = slice2d(xyz[ax2int(our[0])], axis, rad)
            yp = slice2d(xyz[ax2int(our[1])], axis, rad)
            fp = slice2d(f, axis, rad)
            
            N = 200j
            stride = int(N.imag/50)
            cut = cut2-.3
            xcut, ycut = cut, cut
            offset = 1
            mesh_size = lambda ax: 121 if ax == 'y' else 91
            ratio = lambda ax: 1.5 if ax == 'y' else 1
            base_x, base_y = np.linspace(-1, 1, mesh_size(our[0])), np.linspace(-1, 1, mesh_size(our[1]))
            grid_x, grid_y = np.meshgrid(
                np.sign(base_x)*np.abs(base_x)**ratio(our[0])*cut,
                np.sign(base_y)*np.abs(base_y)**ratio(our[1])*cut,
                sparse=False, indexing='ij')
            X, Y = slice1d(xyz[ax2int(our[0])], xyz[ax2int(our[1])], axis, rad)
            #print xyz[0][:,rad,rad], len(xyz[0][:,rad,rad])
            #print xyz[1][rad,:,rad], len(xyz[1][rad,:,rad])
            #print xyz[2][rad,rad,:], len(xyz[2][rad,rad,:])
            fnew = np.nan_to_num(f)
            #print slice2d(fnew, axis, rad)[:,rad], slice2d(fnew, axis, rad)[rad,:]
            spline = RectBivariateSpline(X, Y, slice2d(fnew, axis, rad), kx=3, ky=3, s=smooth)
            grid_z = np.zeros_like(grid_x)
            grid_z.fill(np.NaN)
            mask = grid_x**2 + grid_y**2 < cut**2
            grid_z[mask] = np.maximum(spline(grid_x[mask], grid_y[mask], grid=False), 0)
            zmax = np.round(np.nanmax(grid_z)*1.01, 1-int(np.log(np.nanmax(grid_z))))
            zmin = -1.6*zmax
            print np.nanmin(grid_z), np.nanmax(grid_z), zmax
            ax.set_xlim(-cut2, cut2)
            ax.set_ylim(-cut2, cut2)
            ax.set_zlim(zmin, zmax)
            #surf = ax.plot_wireframe(xp, yp, fp, color='k', linewidth=0.5, rstride=1, cstride=1)
            #surf = ax.plot_wireframe(grid_x, grid_y, grid_z, color='k', linewidth=0.5, rstride=stride, cstride=stride)
            cset = ax.contourf(grid_x, grid_y, grid_z, zdir='z', offset=zmin, cmap=plt.cm.jet,
                levels=np.linspace(0, zmax, 20))
            #cset = ax.contourf(grid_x, grid_y, grid_z, zdir='x', offset=-offset-cut, cmap=plt.cm.jet)
            #cset = ax.contourf(grid_x, grid_y, grid_z, zdir='y', offset=offset+cut, cmap=plt.cm.jet)
            surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=plt.cm.jet, rstride=stride, cstride=stride,
                linewidth=0.25, vmin=0, vmax=zmax, alpha=0.8)
            cbaxes = fig.add_axes([0.82, 0.3, 0.03, 0.6]) 
            cbar = fig.colorbar(surf, cax=cbaxes, drawedges=False)
            cbar.solids.set_edgecolor("face")
            ax.set_xlabel(r'$\xi_%s$' % our[0])
            ax.set_ylabel(r'$\xi_%s$' % our[1])
            break

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)
    axis.set_tick_params(which='both', direction='in', pad=-5)

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

ax.view_init(30, -20)
fig.savefig(outfile)
plt.show()

