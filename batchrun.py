#!/usr/bin/python
# insertparameters.py
# Alex Szatmary
# April 12, 2007
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file

import os, time, os.path

hin = open('batch.parameter','r')

short_tokens = []
tokens = []
list_values = []
n_values = []
m = 0
for line in hin.readlines():
  elements = line.split(',')
  elements[-1] = elements[-1].replace('\n','')
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

files = ['cell', 'fluid', 'force', 'fvs', 'memb', 'rewr', 'visual']
file_in_suffix = '.raw.f'
file_out_suffix = '.run.f'

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
  os.system('cp shpfcta.sph shpfctb.sph shpint.sph sphererad1.out ' + wd)
  fortran_command = 'f90 -o ' + os.path.join(wd, 'cell' + cd)
  for file in files:
      fortran_command = fortran_command + ' ' + os.path.join(wd, file + file_out_suffix)
  os.system(fortran_command)
  os.chdir(wd)
  os.system('nohup ' + 'cell' + cd + ' &')
  os.chdir(os.path.join('..', '..'))
  j = 0
  while i < N_values - 1:
    list_i[j] = list_i[j] + 1
    if list_i[j] == n_values[j]:
      list_i[j] = 0
      j = j + 1
    else:
      break
