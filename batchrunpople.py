#!/usr/bin/python
# insertparameters.py
# Alex Szatmary
# April 12, 2007
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file

import os, time, os.path, sys

hin = open(sys.argv[-1],'r')

short_tokens = []
tokens = []
list_values = []
n_values = []
m = 0
for line in hin.readlines():
  elements = line.split(';')
  elements[-1] = elements[-1].replace('\n','')
  for i in range(len(elements)):
    elements[i] = elements[i].replace('\\n','&\n')
  short_tokens.append(elements[0])
  tokens.append(elements[1])
  list_values.append(elements[2:])
  n_values.append(len(elements)-2)
  m = m + 1

print short_tokens
print tokens
print list_values
print n_values
print m

N_values = 1
for i in n_values:
  N_values = N_values*i

files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
         'math', 'rayTracer', 'meshgen']
file_in_suffix = '.raw.f90'
file_out_suffix = '.run.f90'

visual_files = ['profilemovie']
visual_file_in_suffix = '_raw.m'
visual_file_out_suffix = '_run.m'

list_i = []
values = []
for i in range(m):
  list_i.append(0)
  values.append(0)

for i in range(N_values):
  cd = ''
  for j in range(m):
    values[j] = list_values[j][list_i[j]]
    if short_tokens[j] != '':
      cd = cd+short_tokens[j]+list_values[j][list_i[j]]
  print values
  print cd
  wd = os.path.join('.', 'batch', cd)
  os.mkdir(wd)
  for file in files:
      hin = open(file + file_in_suffix,'r')
      houtcode = open(os.path.join(wd, file + file_out_suffix), 'w')
      for line in hin.readlines():
          for j in range(0,len(tokens)):
            line = line.replace(tokens[j], values[j])
          houtcode.write(line)
      hin.close()
      houtcode.close()
  for file in visual_files:
      hin = open(file + visual_file_in_suffix,'r')
      houtcode = open(os.path.join(wd, file + visual_file_out_suffix), 'w')
      for line in hin.readlines():
          for j in range(0,len(tokens)):
              line = line.replace(tokens[j], values[j])
          houtcode.write(line)
      hin.close()
      houtcode.close()
  hin = open('qsub_pople.raw','r')
  houtcode = open(os.path.join(wd, 'qsub_pople.run'), 'w')
  for line in hin.readlines():
      for j in range(0,len(tokens)):
          line = line.replace(tokens[j], values[j])
      line = line.replace('$cd$', cd)
      houtcode.write(line)
  hin.close()
  houtcode.close()
  print 'compiling'
  fortran_command = 'ifort -o ' + os.path.join(wd, 'cell' + cd)
  for file in files:
      fortran_command = fortran_command + ' ' + os.path.join(wd, file + file_out_suffix)
  fortran_command =  fortran_command + ' -lmpi'
  os.system(fortran_command)
  os.chdir(wd)
  print 'qsubbing'
  os.system('qsub qsub_pople.run')
  print 'qsubbed'
  os.chdir(os.path.join('..', '..'))
  j = 0
  while i < N_values - 1:
    list_i[j] = list_i[j] + 1
    if list_i[j] == n_values[j]:
      list_i[j] = 0
      j = j + 1
    else:
      break
