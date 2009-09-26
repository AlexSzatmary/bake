#!/usr/bin/python
# batchrungfortran.py
# Alex Szatmary
# September 8, 2009
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file. It does 

import os, time, os.path, sys, re


if sys.argv[2] == '-r' or sys.argv[2] == '-rr':
  if sys.argv[2] == '-r':
    task = 'run'
  elif sys.argv[2] == '-rr':
    task = 'rerun'
  system = sys.argv[3]
# Figure out what system I'm running on; make a lot of select cases for this
# Make sure it's a system that has been scripted for, otherwise bad things
# could happen
  if (system != 'gfortran' and system != 'ifort' and 
      system != 'hpc' and system != 'pople'):
    print "Invalid system specified"
    exit(-1)
elif sys.argv[2] == '-e':
  task = 'extract'
  extractfile = sys.argv[3]

#Load bp file
hin = open(sys.argv[1],'r')

short_tokens = []
tokens = []
list_values = []
n_values = []
m = 0
for line in hin.readlines():
  if (line[0] != '#'):
    line = line.replace('\n','').replace('\\n','&\n')
    elements = line.split(';')
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

#Count how many runs I'm going to start
N_values = 1
for i in n_values:
  N_values = N_values*i

#Define which files need tweaking
code_files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
         'math', 'rayTracer', 'meshgen']
file_in_suffix = '.raw.f90'
file_out_suffix = '.run.f90'

#These files are intended to make visualization with Matlab easier.
visual_files = ['profilemovie']
visual_file_in_suffix = '_raw.m'
visual_file_out_suffix = '_run.m'

#m is the number of parameters (not the number of values for the parameters)
#listi and values need to be initialized
list_i = []
values = []
for i in range(m):
  list_i.append(0)
  values.append(0)

pattern = re.compile('\$.+?\$')

tokendict = {}
for i in range(m):
    tokendict[tokens[i]] = i

if task == 'extract':
  hout = open('extract' + extractfile, 'w')
# This is the main loop, setting up each of the runs
for i in range(N_values):
  cd = ''
# Pick the values to be used in this run
  for j in range(m):
    values[j] = list_values[j][list_i[j]]
    if short_tokens[j] != '':
      cd = cd+short_tokens[j]+list_values[j][list_i[j]]
  print values
  print cd

# Do the string replace operations on the values themselves
  for j in range(m):
      foundtoken = re.search(pattern, values[j])
      while foundtoken:
          values[j] = values[j].replace(foundtoken.group(0), 
                                        values[tokendict[foundtoken.group(0)]])
          foundtoken = re.search(pattern, values[j])

  wd = os.path.join('.', 'batch', cd)

  if task == 'run' or task == 'rerun':
    if task == 'run':
      os.mkdir(wd)
  # String replace the tokens for the values
    for file in code_files:
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
    if system == 'hpc':
      hin = open('qsub_script.raw','r')
      houtcode = open(os.path.join(wd, 'qsub_script.run'), 'w')
      for line in hin.readlines():
        for j in range(0,len(tokens)):
          line = line.replace(tokens[j], values[j])
        line = line.replace('$cd$', cd)
        houtcode.write(line)
      hin.close()
      houtcode.close()
    if system == 'pople':
      hin = open('qsub_pople.raw','r')
      houtcode = open(os.path.join(wd, 'qsub_pople.run'), 'w')
      for line in hin.readlines():
        for j in range(0,len(tokens)):
          line = line.replace(tokens[j], values[j])
        line = line.replace('$cd$', cd)
        houtcode.write(line)
      hin.close()
      houtcode.close()

    os.chdir(wd)
    # Insert compile command here  
    if system == 'gfortran':
      fortran_command = 'gfortran -Wall -g -o cell' + cd
    elif system == 'ifort':
      fortran_command = 'ifort -o cell' + cd
    elif system == 'hpc':
      fortran_command = 'mpif90 -o cell' + cd
    elif system == 'pople':
      fortran_command = 'ifort -o cell' + cd
    for file in code_files:
        fortran_command = fortran_command + ' ' + file + file_out_suffix
    os.system(fortran_command)
    if system == 'gfortran' or system == 'ifort':
      os.system('./cell' + cd)
    elif system == 'hpc':
      os.system('qsub qsub_script.run')
    elif system == 'pople':
      os.system('qsub qsub_pople.run')
    os.chdir(os.path.join('..', '..'))

  elif task == 'extract':
    print os.path.join(wd, extractfile)
    if os.path.exists(os.path.join(wd, extractfile)):
      hin = open(os.path.join(wd, extractfile),'r')
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
