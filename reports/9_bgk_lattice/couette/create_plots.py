#!/usr/bin/env python

import os, argparse

parser = argparse.ArgumentParser(description='create a gnuplot script from the template')
parser.add_argument('template', help='template gnuplot script')
parser.add_argument('--name', metavar='str', help='filename without extension')
args, unknown_args = parser.parse_known_args()

for arg in unknown_args:
    assert arg.startswith('--')
    vars(args)[opt] = value

out = '_%s%s' % (args.name, os.path.splitext(args.template)[1])
with open(args.template, 'r') as infile:
    data = infile.read()
    for name, value in args.__dict__.items():
        data = data.replace('<%s>' % name, str(value))
    with open(out, 'w') as outfile:
        outfile.write(data)
    os.chmod(out, 0755)

