#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os, re

#pattern = re.compile('^\$nscale\$')

for root, dirs, files in os.walk('./bp/mix'):
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
            # prolate
            if 'dsqrt(0.5d0)' in line:
                line = line.replace('dsqrt(0.5d0)', 
                                    '0.65368179876683d0')
            # oblate
            if 'dsqrt(0.75d0)' in line:
                line = line.replace('dsqrt(0.75d0)', 
                                    '0.83071445098496d0')
            hout.write(line)
        hout.close()
