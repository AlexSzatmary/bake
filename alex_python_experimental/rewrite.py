#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os

killlines = ['# ellipa is the semi-major axis length in the n1 direction\n',
             '# b in the n2, and c in the n3\n',
             '# Each ellipsoid has three entries in each of these arrays,\n',
             '# ie, entries 1:3 in $ellipn2$ correspond to the n2 vector\n',
             '# of the first ellipse, entries 4:6 correspond to the n2\n',
             '# vector of the second ellipse, and so on.\n',
             '$ellipn1$;(/0/)\n',
             '$ellipn2$;(/0/)\n',
             '$ellipn3$;(/0/)\n',]

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
            if '# Specify geometry for ellipses' in line:
                hout.write(
                    """# Specify geometry for ellipses
# ellipa, b, c describe the semi-principal axes; each is a vector going from 
# the center to semi-principal points. These vectors should be orthogonal.
# Each ellipsoid has three entries in each of these arrays,
# ie, entries 1:3 in $ellipb$ correspond to the b vector
# of the first ellipse, entries 4:6 correspond to the b
# vector of the second ellipse, and so on.
""")
            else:
                flag = 0
                for killline in killlines:
                    if killline == line:
                        flag = 1
                        break
                if flag == 0:
                    hout.write(line)
        hout.close()
