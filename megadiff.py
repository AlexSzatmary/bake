#!/usr/bin/python
# megadiff.py
# Alex Szatmary
# 2009-08-08
# This is useful in comparing different runs with the same set of
# parameters. The idea is that, while revising the code, megadiff can
# be used when numerically identical results are expected from run to
# run. This is useful when making stylistic changes to the code, but
# not when implementing new physics.
# To use, copy a directory from the batch directory to the main
# directory of the repository, so you have a tree that looks something
# like this:
# cell/batch/myjob
# cell/myjob
# Run this script:
# ./megadiff.py myjob
# This will then diff a handful of files. If you see noise, you broke
# something affecting the numerics.

import os, time, os.path, sys

wd = sys.argv[-1]

files = ['TaylorDF__00001.txt', 'capsulev__00001.txt', 
         'capsulex__00001.txt', 'checkinit.txt', 'meanfluidv.txt', 'meanforce.txt',
         'cappro0001_000000.txt', 'solidforce00000.txt', 
         'solidforce00005.txt',
         'solidforce00010.txt', 'solidnodes00000.txt', 
         'solidnodes00005.txt', 'solidnodes00010.txt', 'status.txt', 
         'uvwpdump__00000.txt', 'uvwpdump__00005.txt',
         'uvwpdump__00010.txt', 'wprofile__00000.txt', 
         'wprofile__00001.txt', 'wprofile__00002.txt', 
         'wprofile__00003.txt', 'wprofile__00004.txt',
         'wprofile__00005.txt', 'wprofile__00006.txt',
         'wprofile__00007.txt', 'wprofile__00008.txt',
         'wprofile__00009.txt', 'wprofile__00010.txt',
         'minmaxx___00001.txt', 'minmaxy___00001.txt', 'minmaxz___00001.txt',
         'volumearea00001.txt',
         'thumbprint.txt']
for file in files:
    os.system('diff batch/'+wd+'/'+file+' check/'+wd+'/'+file)
