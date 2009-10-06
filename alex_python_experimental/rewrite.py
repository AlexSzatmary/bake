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
        linesout = []
        label = ''
        for line in lines:
            if line[0] != '#':
                l = line.split(';')
                if l[0] != '':
                    label += l[0]+l[1]
                linesout.append(';'.join(l[1:]))
            else:
                linesout.append(line)
        hout.write('$label$;'+label+'\n')
        for line in linesout:
            hout.write(line)
        hout.close()
