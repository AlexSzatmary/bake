#!/usr/bin/python
# batchrun.py
# Alex Szatmary
# September 8, 2009
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file. It does 

import os, time, os.path, sys, re, listruns

i = 0

while i < len(sys.argv):
  try:
    if '.py' in sys.argv[i]:
      pass
    elif sys.argv[i] == '-r' or sys.argv[i] == '-rr':
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      if sys.argv[i] == '-r':
        task = 'run'
      elif sys.argv[i] == '-rr':
        task = 'rerun'
      i += 1
      system = sys.argv[i]
# Figure out what system I'm running on; make a lot of select cases for this
# Make sure it's a system that has been scripted for, otherwise bad things
# could happen
      if (system != 'gfortran' and system != 'gfortranopenmp' and 
          system != 'ifort' and system != 'hpc' and system != 'pople'):
        raise Exception('Invalid system specified')
    elif sys.argv[i] == '-e':
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      task = 'extract'
      i += 1
      if '-sn' == sys.argv[i]:
        shortname = 1
        i += 1
      else:
        shortname = 0
      extractfile = sys.argv[i]
      if shortname == 1:
        if extractfile == 'DF':
          extractfile = 'TaylorDF__00001.txt'
        else:
          raise Exception('Invalid shortname')
# Perform operation on a Slice of the runs      
    elif sys.argv[i] == '-s':
      if 'slice_start' in dir() or 'slice_end' in dir():
        raise Exception('Multiple slices specified')
      i += 1
      if '-' in sys.argv[i]:
        (slice_start, slice_end) = sys.argv[i].split('-')
      else:
        slice_start = 0
        slice_end = int(sys.argv[i])
      if slice_start == '':
        slice_start = 0
      if slice_end == '':
        slice_end = 0
      slice_start = int(slice_start)
      slice_end = int(slice_end)
    elif sys.argv[i] == '-l':
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      task = 'list'
    else:
      if 'myfile' in dir():
        raise Exception('Batch parameter file already specified')
      myfile = sys.argv[i]
    i += 1
  except Exception, data:
    if data[0] == 'Invalid system specified':
      print data[0]
      exit(-1)
    elif data[0] == 'Multiple tasks requested':
      print data[0]
      exit(-2)
    elif data[0] == 'Batch parameter file already specified':
      print data[0]
      exit(-3)
    elif data[0] == 'Invalid shortname':
      print data[0]
      exit(-4)
    else:
      print 'I don\'t know what exception I\'m dealing with.'
      raise
      exit(-100)


###############################################################################
# End processing of command line parameters
###############################################################################

(tokens, list_values, n_values, N_values, m) = listruns.LoadBPFile(myfile)

print tokens
print list_values
print n_values
print m
print 'Number of runs in file ', N_values

if 'slice_start' not in dir():
  slice_start = 0

if 'slice_end' not in dir() or slice_end == 0:
  slice_end = N_values

if task == 'list':
  exit(0)

#Define which files need tweaking
code_files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
         'math', 'rayTracer', 'meshgen']
file_in_suffix = '.raw.f90'
file_out_suffix = '.run.f90'

#These files are intended to make visualization with Matlab easier.
visual_files = ['profilemovie']
visual_file_in_suffix = '_raw.m'
visual_file_out_suffix = '_run.m'

pattern = re.compile('\$.+?\$')

tokendict = {}
for i in range(m):
    tokendict[tokens[i]] = i


if task == 'extract':
  hout = open('extract' + extractfile, 'w')

# This is the main loop, setting up each of the runs
for values in listruns.ItRunValues(list_values, tokens, n_values, N_values, m, 
                              pattern, tokendict, slice_start, slice_end):
# Do the string replace operations on the values themselves
  cd = values[tokendict['$label$']]
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
    elif system == 'gfortranopenmp':
      fortran_command = 'gfortran -Wall -fopenmp -g -o cell' + cd
    elif system == 'ifort':
      fortran_command = 'ifort -o cell' + cd
    elif system == 'hpc':
      fortran_command = 'mpif90 -o cell' + cd
    elif system == 'pople':
      fortran_command = 'ifort -o cell' + cd
    for file in code_files:
        fortran_command = fortran_command + ' ' + file + file_out_suffix
    os.system(fortran_command)
    if system == 'gfortran' or system == 'gfortranopenmp' or system == 'ifort':
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

if task == 'extract':
  hout.close()
  hin = open('extract' + extractfile, 'r')
  print hin.read()
  hin.close()
