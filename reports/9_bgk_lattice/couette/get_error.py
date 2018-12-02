#!/usr/bin/env python
import sys, argparse
import numpy as np

parser = argparse.ArgumentParser(description='Calculate deviation from the reference value')
parser.add_argument('txtfile', type=str, help='filename with data')
parser.add_argument('-c', '--column', type=int, default=2, help='the column in file')
parser.add_argument('-v', '--value', type=float, default=0, help='the reference value')
parser.add_argument('-U', type=float, default=0.02, help='the speed in the Couette-flow problem')
args = parser.parse_args()

data = np.loadtxt(args.txtfile)
print('{:.7f}'.format(data[0][args.column]/args.U - args.value))

