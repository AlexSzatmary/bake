#!/usr/bin/python

#import optparse
#optparser = optparse.OptionParser()
#options, arguments = optparser.parse_args()
#print arguments
#for argument in arguments:
#    print argument

hinstretches = open('stretches_01000.txt', 'r')
hinshpa = open('shpa0001.txt', 'r')
W = 0.
area = 0.
for line in hinstretches:
    elements = line.split(' ')
    while '' in elements:
        elements.remove('')
    l1 = float(elements[1])
    l2 = float(elements[2][:-1])
    line = hinshpa.readline()
    elements = line.split(' ')
    while '' in elements:
        elements.remove('')
    area_i = float(elements[3])
    area = area + area_i
    W = W + (l1*l1 + l2*l2 + 1/(l1*l1*l2*l2) - 3)*area_i/3
W = W/area
hinstretches.close()
hinshpa.close()
print '%.14e'%W
