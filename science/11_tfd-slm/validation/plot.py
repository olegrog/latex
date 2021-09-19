#!/usr/bin/env python

import matplotlib, argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='plot with error bars for a set of powers')
parser.add_argument('pdffile', help='an output PDF file')
parser.add_argument('-x', default='power/speed', help='value for the x axis')
parser.add_argument('-y', default='width', help='value for the y axis')
parser.add_argument('--xlabel', default='P/U, J/mm', help='label for the x axis')
parser.add_argument('--ylabel', default='d, µm', help='label for the y axis')
parser.add_argument('--logx', action='store_true', help='log scaling on the x axis')
parser.add_argument('--logy', action='store_true', help='log scaling on the y axis')
parser.add_argument('--use-tex', action='store_true', help='use tex for text')
args = parser.parse_args()

matplotlib.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.linewidth': 1,
    'errorbar.capsize': 3,
    'text.usetex': args.use_tex,
    'pgf.rcfonts': False,
    'font.family': 'serif',
    'pgf.preamble': '\n'.join([
         r'\usepackage{physics, siunitx}',
    ]),
    'figure.figsize': [6, 4]
})
if args.use_tex:
    matplotlib.use('pgf')

def set_colors():
    # https://stackoverflow.com/a/55652330/2531400
    cmap = matplotlib.cm.get_cmap('Set1')
    axes = plt.gca()
    axes.set_prop_cycle(color=cmap.colors)

def load_data(filename, columns):
    data = np.loadtxt(filename).T
    result = dict(zip(columns, data))
    result['power/speed'] = result['power']/result['speed']
    return result

def analyze_data(data):
    result, err = {}, {}
    for p in np.unique(data['power']):
        for s in np.unique(data['speed']):
            mask = np.argwhere((data['power'] == p) & (data['speed'] == s))
            if not np.sum(mask):
                continue
            for key, value in data.items():
                if key not in result:
                    result[key], err[key] = [], []
                result[key].append(np.mean(value[mask]))
                err[key].append(np.std(value[mask]))
    for key in result:
        result[key] = np.array(result[key])
        err[key] = np.array(err[key])
    return result, err

def plot_profile(data, fmt, label, err=None, lw=1, ms=3, **kwargs):
    set_colors()
    if args.y not in data:
        return
    x, y = data[args.x], data[args.y]
    for p in np.unique(data['power']):
        mask = np.argwhere(data['power'] == p)[:,0]
        yerr = err if err is None else err[args.y][mask]
        plt.errorbar(x[mask], y[mask], yerr=yerr, fmt=fmt, lw=lw, markersize=ms,
            label=label + f' ({p:.0f}W)', **kwargs)

# Experimental data 1
data = load_data('exper1.txt', ['power', 'speed', 'width'])
plot_profile(data, '^', 'Experiment1')

# Experimental data 2
data = load_data('exper2.txt', ['power', 'speed', 'depth', 'width', 'height'])
#plot_profile(data, 'v', 'Experiment2')
data, err = analyze_data(data)
plot_profile(data, '-', 'Experiment2', err=err)

# Simulation data
data = load_data('numer1.txt', ['power', 'speed', 'width', 'depth'])
plot_profile(data, 'o', 'Simulation1', ms=5, mfc='none')

#plt.xlim(np.min(x), np.max(x))
#plt.ylim(np.min(y), np.max(y))

plt.xlabel(args.xlabel)
plt.ylabel(args.ylabel)
plt.legend(loc="upper left")
P, U = 113, 700 # Recommended for SS316L by Trumpf
plt.axvline(x=P/U, c='b', ls=':', lw=.5)

if args.logx:
    plt.semilogx()
if args.logy:
    plt.semilogy()

plt.savefig(args.pdffile, bbox_inches='tight', transparent=True)
