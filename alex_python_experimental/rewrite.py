#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os

for root, dirs, files in os.walk('./bp/'):
    if '.svn' in dirs:
	dirs.remove('.svn')
    for fname in files:
	if '~' in fname:
	    files.remove(fname)
    if 'parameters.txt' in files:
        files.remove('parameters.txt')
    print root, dirs, files
    for fname in files:
        hin = open(os.path.join(root, fname))
        lines = hin.readlines()
        hin.close()
        hout = open(os.path.join(root, fname), 'w')
        for line in lines:
            if '$ncap$' not in line:
                hout.write(line)
            else:
                hout.write(line.replace('$ncap$', '$nsph$'))
                hout.write('$ncap$;($nsph$+$nellip$)\n')
                hout.write('$nellip$;0\n')
        hout.write('\n# Specify geometry for ellipses\n')
        hout.write('# ellipa is the semi-major axis length in the n1 direction\n')
        hout.write('# b in the n2, and c in the n3\n')
        hout.write('$ellipa$;(//)\n')
        hout.write('$ellipb$;(//)\n')
        hout.write('$ellipc$;(//)\n')
        hout.write('# Each ellipsoid has three entries in each of these arrays,\n')
        hout.write('# ie, entries 1:3 in $ellipn2$ correspond to the n2 vector\n')
        hout.write('# of the first ellipse, entries 4:6 correspond to the n2\n')
        hout.write('# vector of the second ellipse, and so on.\n')
        hout.write('$ellipn1$;(//)\n')
        hout.write('$ellipn2$;(//)\n')
        hout.write('$ellipn3$;(//)\n')
        hout.close()
