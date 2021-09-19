#!/usr/bin/python

import numpy as np
import sys

b, a = 2, 1
err = 1e-5

def term(x, i):
    t, tt = 1 + ((b-a)*x-x*x/2) / (a*b), ((b-a)-x) / (a*b)
    alpha, der = np.arccosh(t), tt/np.sqrt(t*t-1)
    A = b*np.sinh((i+1)*alpha) - a*np.sinh(i*alpha)
    AA = b*(i+1)*np.cosh((i+1)*alpha) - a*i*np.cosh(i*alpha)
    return -der / A * (np.cosh(alpha) - np.sinh(alpha) * AA / A)

def func(x):
    result = np.zeros(x.shape)
    for i in range(1, 10**4):
        t = term(x, i)
        result += t
        if np.any(result == 0):
            continue
        if np.all(np.abs(t/result) < err):
            break
    return result, i

X = np.logspace(np.log10(1e-3), np.log10(0.95), num=30)
Y, i = func(X)
print("num of iters = %d" % i)
    
np.savetxt(sys.argv[1], np.transpose((X, Y**-1)), fmt='%1.4e') 
