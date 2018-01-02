#!/usr/bin/env python
    
import glob, os, argparse
import pandas as pd
import pylab as py

parser = argparse.ArgumentParser(description='CPU analyser for the plane Couette-flow problem')
parser.add_argument('dirs', nargs='+', help='list of directories')
args = parser.parse_args()

def read_log(filename):
    N, m1, m2 = os.path.basename(filename).split('.')[0].split('-')
    result = { 'N': N, 'model1': m1, 'model2': m2 }
    with open(filename, 'r') as lines:
        for line in lines:
            if 'User time' in line:
                result['cpu'] = float(line.split()[-1])
    return result

def update(dst, src, key):
    if key in dst:
        dst[key] += src
    else:
        dst[key] = src

def scheme_name(result):
    if result['model1'] != result['model2']:
        return 'hybrid'
    if result['model1'] == 'dvm':
        return 'dvm'
    else:
        return 'lbm'

data = {}
for dirname in args.dirs:
    for filename in glob.glob(os.path.join(dirname, '*.log')):
        print filename
        result = read_log(filename)
        print result
        N = int(result['N'])
        if not N in data:
            data[N] = {}
        update(data[N], result['cpu'], scheme_name(result))
        
df = pd.DataFrame(data).transpose()
py.plot(df.dvm/df.hybrid, '*-')
py.loglog()
py.xlabel('Number of cells in the Knudsen layer')
py.ylabel('Acceleration (t_DVM/t_hybrid)')
py.grid()
py.xlim(2, 1e3)
py.ylim(1, 15)
py.show()

