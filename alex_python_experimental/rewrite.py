#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os

for root, dirs, files in os.walk('./bp/test'):
    if '.svn' in dirs:
	dirs.remove('.svn')
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
            if ';$nstep$;10' in line:
                line = line.replace(';$nstep$;10', ';$nstep$;5')
            hout.write(line)
        hout.close()
