#!/usr/bin/python

# This reads a solution.txt file exported by poisson.py and compares that
# solution to an ideal one given in a bake parameter file. The L_infinity norm
# is output; this is the maximum difference between the computed and ideal
# solutions. solution.txt must already exist before running this script.

from model import *
import math


# Gets an ideal solution function from bake
def ideal(x):
    return @ideal@

hin = open('solution.txt')

table = []

for l in hin:
    (x, u) = [float(s) for s in l.split()]
    table.append([x, u, ideal(x), abs(ideal(x) - u)])
#    print('%14.7e %14.7e %14.7e %14.7e' % (x, u, ideal(x), abs(ideal(x) - u)))

print max(zip(*table)[3])

hin.close()
