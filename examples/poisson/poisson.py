#!/usr/bin/python

# kludgey way to get bake parameters into code
from model import *
import math

def tridiag_solver(b):
    """Solves the tridiagonal system diag(1,-2,1)x = b, given b"""
    for i in xrange(1, len(b)):
        b[i] = b[i] + b[i - 1] * i/(i + 1.)
    x = [0. for i in xrange(len(b))]
    x[-1] = - b[-1] * (i + 1.)/(i + 2.)
    for i in xrange(len(b) - 2, -1, -1):
        x[i] = -(b[i] - x[i + 1])*(i + 1)/(i + 2)
    return x

def set_up():
    # Set up b in terms of its forcing function
    b = [f(x_l + h*i)*h*h for i in xrange(1, n - 1)]
    # Apply boundary conditions to b
    b[0] -= u_l
    b[-1] -= u_r

    u = tridiag_solver(b)

    hout = open('solution.txt', 'w')
    hout.write('%21.14e %21.14e\n' % (x_l, u_l))
    for i in xrange(0, n - 2):
        hout.write('%21.14e %21.14e\n' % (x_l + (i + 1)*h, u[i]))
    hout.write('%21.14e %21.14e\n' % (x_r, u_r))

if __name__== "__main__":
    set_up()
