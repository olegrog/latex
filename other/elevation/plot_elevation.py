#!/usr/bin/env python3

import sys, argparse, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import gpxpy
import gpxpy.gpx

str2pair = lambda s: [float(item) for item in s.split(':')]

parser = argparse.ArgumentParser(description='Plat an elevation profile')
parser.add_argument('input', type=str, help='input filename')
parser.add_argument('-o', '--output', type=str, default='elevation.pdf', help='output filename')
parser.add_argument('-m', '--markers', action='store_true', help='add markers at the peaks')
parser.add_argument('-s', '--figsize', type=str2pair, default='8:3', help='figure size')
parser.add_argument('-t', '--tick-spacing', type=float, default=10, help='tick spacing for x-axis')
parser.add_argument('-w', '--width', type=float, default=10, help='minumum width of peaks (points)')
parser.add_argument('-p', '--prominence', type=float, default=100, help='minumum prominence of peaks (m)')
parser.add_argument('--pdf', action='store_true', help='save a PDF file instead')
args = parser.parse_args()

class Style:
    line = { 'linestyle': '-', 'linewidth': 1, 'color': 'k' }
    dashed = { 'linestyle': '--', 'linewidth': 0.5, 'color': 'k' }
    point = { 'color': 'black', 'marker': '.', 'linestyle': 'None', 'zorder': 10 }
    annotate = { 'xytext': (0,5), 'textcoords': 'offset points', 'ha': 'center' }

params = {
    'axes.labelsize': 11,
    'font.size': 11,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.figsize': args.figsize
}
plt.rcParams.update(params)

if args.pdf:
    matplotlib.use('pgf')
    plt.rcParams['text.usetex'] = True
    plt.rcParams['pgf.rcfonts'] = False

with open(args.input, 'r') as gpx_file:
    gpx = gpxpy.parse(gpx_file)

X, Y = [], []
X_days = []
previous_point = None

for track in gpx.tracks:
    for segment in track.segments:
        for point in segment.walk(only_points=True):
            X.append(X[-1] + point.distance_2d(previous_point)/1e3 if previous_point else 0)
            Y.append(point.elevation);
            previous_point = point
    X_days.append(X[-1])
print(f'Total length = {X[-1]:.1f} km')

# 1. Draw curves
plt.plot(X, Y, **Style.line)
plt.fill_between(X, Y, np.min(Y), color='lightgray')
for x in X_days[:-1]:
    plt.axvline(x, **Style.dashed)

# 2. Annotate peaks
X, Y = np.array(X), np.array(Y)
peaks, _ = find_peaks(Y, width=args.width, prominence=args.prominence)
print(f'Total peaks found = {len(peaks)}')
if args.markers:
    plt.plot(X[peaks], Y[peaks], **Style.point)
for i, peak in enumerate(peaks):
    print(f'{i+1}) L = {X[peak]:.1f} km, H = {Y[peak]} m')
    plt.annotate(f'{i+1}', (X[peak], Y[peak]), **Style.annotate)

# 3. Adjust margins and ticks
plt.autoscale(enable=True, axis='x', tight=True)
delta = 0.1*(np.max(Y) - np.min(Y))
plt.ylim(np.min(Y), np.max(Y) + delta)
ax = plt.gca()
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(args.tick_spacing))

filename = args.output if args.output else f'{args.mode}.pdf'
plt.tight_layout(pad=1)
if args.pdf:
    plt.savefig(filename, bbox_inches='tight')
else:
    plt.show()

