#!/usr/bin/python

import optparse

optparser = optparse.OptionParser()

options, arguments = optparser.parse_args()

print arguments

for argument in arguments:
    print argument
    hin = open(argument, 'r')
    DFlist = []
    for line in hin.readlines():
	elements = line.split(' ')
	while '' in elements:
	    elements.remove('')
#	print ' '.join(elements)[:-1]
	DFlist.append(float(elements[4]))
    hin.close()
#    print DFlist
    n = len(DFlist)

    if 'sp' not in argument:
	print 'tweaking DF'
	for i in xrange(n):
	    DFlist[i] = (3*DFlist[i]-1)/(3-DFlist[i])

    DFdifflist = []
    for i in xrange(n):
	DFdifflist.append(abs((DFlist[i] - (DFlist[-1]*i/(n-1) + 
		   DFlist[0]*(n-i-1)/(n-1)))/DFlist[i]))

#    for diff in DFdifflist:
#	print diff
    print max(DFdifflist)
