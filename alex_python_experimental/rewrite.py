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
            if '$a_prestress$,' not in line:
                hout.write(line)
            else:
                hout.write(line.replace('$a_prestress$,', '$a_prestress$;'))
        hout.close()
