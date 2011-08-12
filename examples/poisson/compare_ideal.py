#!/usr/bin/python

from model import *
import math

def ideal(x):
    return @ideal@

hin = open('solution.txt')

table = []

for l in hin:
    (x, u) = [float(s) for s in l.split()]
#    print('%14.7e %14.7e %14.7e %14.7e' % (x, u, ideal(x), abs(ideal(x) - u)))
    table.append([x, u, ideal(x), abs(ideal(x) - u)])

print max(zip(*table)[3])

hin.close()
