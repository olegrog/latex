#!/usr/bin/env python

import os, argparse

parser = argparse.ArgumentParser(description='create a gnuplot script from the template')
parser.add_argument('template', help='template gnuplot script')
parser.add_argument('--name', metavar='str', help='filename without extension')
parser.add_argument('--xlabel', metavar='str', help='x label')
parser.add_argument('--ylabel', metavar='str', help='y label')
parser.add_argument('--patch', metavar='str', help='patch name')
parser.add_argument('--column', type=int, metavar='value', help='column to plot')
parser.add_argument('--xcoord', type=float, default=1, metavar='value', help='x-coord of the y-label')
parser.add_argument('--ytics', type=float, default=1, metavar='value', help='step of ytics')
parser.add_argument('--kn', type=float, default=1, metavar='value', help='a Knudsen number for transform (x,y)-coord')
args = parser.parse_args()

out = '_%s%s' % (args.name, os.path.splitext(args.template)[1])
with open(args.template, 'r') as infile:
    data = infile.read()
    for name, value in args.__dict__.items():
        data = data.replace('<%s>' % name, str(value))
    with open(out, 'w') as outfile:
        outfile.write(data)
    os.chmod(out, 0755)

