#!/usr/bin/python
# extractDF.py
# Alex Szatmary
# March 27, 2008
# Extracts DF data from all runs performed.

import os, time, os.path, sys

hin = open(sys.argv[-1],'r')

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

files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual']
file_in_suffix = '.raw.f'
file_out_suffix = '.run.f'

list_i = []
values = []
for i in range(m):
  list_i.append(0)
  values.append(0)

hout = open('extractDF.txt', 'w')
for i in range(N_values):
  cd = ''
  for j in range(m):
    values[j] = list_values[j][list_i[j]]
    if short_tokens[j] != '':
      cd = cd+short_tokens[j]+list_values[j][list_i[j]]
  wd = os.path.join('.', 'batch', cd)
  print os.path.join(wd, 'fort.204')
  if os.path.exists(os.path.join(wd, 'fort.204')):
    hin = open(os.path.join(wd, 'fort.204'),'r')
    hout.write(wd + ',' + hin.readlines()[-1])
    hin.close()
  j = 0
  while i < N_values - 1:
    list_i[j] = list_i[j] + 1
    if list_i[j] == n_values[j]:
      list_i[j] = 0
      j = j + 1
    else:
      break
