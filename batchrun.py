#!/usr/bin/python
# batchrun.py
# Alex Szatmary
# September 8, 2009
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file. It does 

import os, time, os.path, sys, re, listruns, optparse

optparser = optparse.OptionParser()
optparser.add_option('--plot', '-p')
optparser.add_option('--run', '-r')
optparser.add_option('--rerun', '-R')
optparser.add_option('--extract', '-e')
optparser.add_option('--slice', '-s')
optparser.add_option('--list', '-l', action='store_true')
optparser.add_option('--overwrite', '-o', action='append',
		     help="Overwrite a line in a batch parameter file")
optparser.add_option('--foreach', '-E')
options, arguments = optparser.parse_args()

print options
print arguments

try:
  if options.extract:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'extract'
    extractfile = options.extract
    if extractfile == 'DF':
      extractfile = 'TaylorDF__00001.txt'

  if options.run or options.rerun:
    if 'task' in dir() or (options.run and options.rerun):
      raise Exception('Multiple tasks requested')
    if options.run:
      task = 'run'
      system = options.run
    elif options.rerun:
      task = 'rerun'
      system = options.rerun
    # Figure out what system I'm running on; make a lot of select cases for this
    # Make sure it's a system that has been scripted for, otherwise bad things
    # could happen
    if (system != 'gfortran' and system != 'gfortranopenmp' and 
        system != 'ifort' and system != 'hpc' and system != 'pople'
        and system != 'sun'):
      raise Exception('Invalid system specified')

  # Perform operation on a Slice of the runs      
  if options.slice:
    if '-' in options.slice:
      (slice_start, slice_end) = options.slice.split('-')
    else:
      slice_start = 0
      slice_end = int(options.slice)
    if slice_start == '':
      slice_start = 0
    if slice_end == '':
      slice_end = 0
    slice_start = int(slice_start)
    slice_end = int(slice_end)

  if options.list:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'list'

  if options.plot:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'plot'
    extractfile = options.plot
    if extractfile == 'DF':
      extractfile = 'TaylorDF__00001.txt'

  if options.foreach:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'foreach'
  myfile = arguments[0]
  print task
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

if options.overwrite:
  lines = options.overwrite
else:
  lines = []

print lines
hin = open(myfile,'r')
lines += hin.readlines()
hin.close()

(tokens, list_values, n_values, N_values, tokendict, m) = \
    listruns.parseBPlines(lines)

print tokens
print list_values
print n_values
print m
print 'Number of runs in file ', N_values

if 'slice_start' not in dir():
  slice_start = 0

if 'slice_end' not in dir() or slice_end == 0:
  slice_end = N_values


#Define which files need tweaking
code_files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
         'math', 'meshgen']
file_in_suffix = '.raw.f90'
file_out_suffix = '.run.f90'

#These files are intended to make visualization with Matlab easier.
visual_files = ['profilemovie']
visual_file_in_suffix = '_raw.m'
visual_file_out_suffix = '_run.m'

pattern = re.compile('\$.+?\$')

# tokendict = {}
# for i in range(m):
#     tokendict[tokens[i]] = i


if task == 'extract' or task == 'plot':
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
    elif system == 'sun':
      fortran_command = 'f95 -o cell' + cd
    elif system == 'ifort':
      fortran_command = 'ifort -o cell' + cd
    elif system == 'hpc':
      fortran_command = 'mpif90 -o cell' + cd
    elif system == 'pople':
      fortran_command = 'ifort -o cell' + cd
    for file in code_files:
        fortran_command = fortran_command + ' ' + file + file_out_suffix
    os.system(fortran_command)
    if system == 'gfortran' or system == 'gfortranopenmp' or system == 'ifort'\
           or system == 'sun':
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
  elif task == 'plot':
    print os.path.join(wd, extractfile)
    if os.path.exists(os.path.join(wd, extractfile)):
      hin = open(os.path.join(wd, extractfile),'r')
      mix = values[tokendict['$mix$']]
      capillary_no = values[tokendict['$capillary_no$']]
      line = mix + ' ' + capillary_no + hin.readlines()[-1]
      hout.write(wd + ' ' + line)
      hin.close()
  elif task == 'list':
    print cd
  elif task == 'foreach':
    os.chdir(wd)
    os.system(options.foreach)
    os.chdir(os.path.join('..', '..'))

if task == 'extract' or task == 'plot':
  hout.close()
  hin = open('extract' + extractfile, 'r')
  print hin.read()
  hin.close()
