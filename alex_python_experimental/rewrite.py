#!/usr/bin/python
# rewrite.py
# I made this Python script to make it easy to batch rewrite my bp files

import os

for root, dirs, files in os.walk('./bp'):
    print dirs
    if '.svn' in dirs:
	dirs.remove('.svn')
    for fname in files:
	if '~' in fname:
	    files.remove(fname)
    print root, dirs, files
