#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os

for root, dirs, files in os.walk('./bp/'):
    if '.svn' in dirs:
	dirs.remove('.svn')
    if 'test' in dirs:
        dirs.remove('test')
    for fname in files:
	if '~' in fname:
	    files.remove(fname)
    print root, dirs, files
    for fname in files:
        hin = open(os.path.join(root, fname))
        lines = hin.readlines()
        hin.close()
        hout = open(os.path.join(root, fname), 'w')
        for line in lines:
            hout.write(line)
            if ';$nstep$;' in line:
                line = line.replace(';$nstep$;', ';$nstep$;')
                hout.write(line)
        hout.close()
