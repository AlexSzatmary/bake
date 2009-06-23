#!/usr/bin/python
# insertparameters.py
# Alex Szatmary
# April 12, 2007
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file

import os, time

tokens = ['$nsnode$', '$nselm$', '$lngx$', '$Eh$', '$capillary_no$', '$FVS$']
values = ['10242', '20480', '5', '1.', '1.d-2', '1']
files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual']
file_in_suffix = '.raw.f'
file_out_suffix = '.run.f'

for file in files:
    hin = open(file + file_in_suffix,'r')
    houtcode = open(file + file_out_suffix, 'w')
    for line in hin.readlines():
        for j in range(0,len(tokens)):
            line = line.replace(tokens[j], values[j])
        houtcode.write(line)
    hin.close()
    houtcode.close()
