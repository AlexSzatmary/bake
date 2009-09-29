#!/usr/bin/python

import sys
hin = open(sys.argv[1], 'r')
lines = hin.readlines()
hin.close()

x = []
y = []
z = []

for line in lines:
    elm = line.replace('\n', '').split(' ')
    while '' in elm:
	elm.remove('')
    x.append(float(elm[0]))
    y.append(float(elm[1]))
    z.append(float(elm[2]))

print min(x), ' < x < ', max(x)
print min(z), ' < y < ', max(y)
print min(z), ' < z < ', max(z)
