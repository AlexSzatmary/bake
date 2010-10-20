#!/usr/bin/python
# batchrun.py
# Alex Szatmary
# September 8, 2009
# This script takes a raw Fortran file (By raw, I mean not ready to compile,
# having token codes like '$lngx$' present instead of an actual value, like 6.)
# and translates each token to a value set in this file. It does 

import os, os.path, sys, re, listruns, optparse

optparser = optparse.OptionParser()
optparser.add_option('--plot', '-p',
                     help="Makes a datafile for a DF-m plot. "
                     "This is only useful right now for Alex's research, but "
                     "if you want code for something like this, let him know.")
optparser.add_option('--run', '-r',
                     help="""Start a run; specify a system to run on.
                     eg, "-r foo" asks to run the code on system foo.
                     A system must be specified. Currently used systems
                     include "hpc", "pople", "tara", and "sun".""")
optparser.add_option('--rerun', '-R',
                     help="""Like -r, but restarts a stopped run.
                     The only difference between this and -r is that new
                     directories are not made.""",)
optparser.add_option('--extract', '-e',
                     help="""Extracts last line of the data file specified
                     here, for each run, prints these lines, and saves them.
                     """)
optparser.add_option('--slice', '-s',
                     help="""Selects a subset of the runs specified in the
                     file, eg, -s 5:9 does runs 5, 6, 7, and 8 out of however
                     many runs would be referred to in the given file.""")
optparser.add_option('--list', '-l', action='store_true',
                     help="""Lists the jobs that would be operated on with the
                     given parameter file and options.""")
optparser.add_option('--overwrite', '-o', action='append',
		     help="""Overwrite a line in a batch parameter file,
                     eg, "-o '\$foo\$;bar;baz'" replaces a parameter line
                     starting with "$foo$" with "$foo$;bar;baz".
                     This option can be used repeatedly.
                     (Note: if the parameter specified is absent from the file,
                     the new line will simply be added to the options in the
                     file, it won't overwrite anything.)
                     """)
optparser.add_option('--foreach', '-E',
                     help="""Execute a command in each job specified, eg, "tail
                     TaylorDF__0001.txt"
                     """)
optparser.add_option('--fit', '-f', action='store_true',
                     help="""Does a curve fit to 1-exp(t) using gnuplot for
                     Taylor DF.""")
optparser.add_option('--clone', '-c', 
                     help="""Copies one file from each specified directory to
                     a clone in Alex/clone.
                     """)
optparser.add_option('--backup', '-b', action='store_true',
                     help="""Backs up selected runs to Alex/backup. This is
                     helpful to do between runs, after checking to make sure
                     that the checkpoint and data files are saved safely.
                     """)
optparser.add_option('--restore', '-t', action='store_true',
                     help="""Brings back runs backed up by the backup option.
                     """)
optparser.add_option('--timeseries', '-i', action='store_true',
                     help="""Make time-series DF plots for selected runs
                     """)
options, arguments = optparser.parse_args()

#print options
#print arguments

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
    system_list = ['gfortran', 'gfortranopenmp', 'ifort' ,'hpc', 'pople',
                   'tara', 'sun', 'popledebug']
    if system not in system_list:
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

  if options.fit:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'fit'

  if options.foreach:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'foreach'

  if options.clone:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'clone'
    clonefile = options.clone

  if options.backup:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'backup'

  if options.restore:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'restore'

  if options.timeseries:
    if 'task' in dir():
      raise Exception('Multiple tasks requested')
    task = 'timeseries'

  myfile = arguments[0]
#  print task

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

# print lines
hin = open(myfile,'r')
lines += hin.readlines()
hin.close()

(tokens, list_values, n_values, N_values, tokendict, m) = \
    listruns.parseBPlines(lines)

#print tokens
#print list_values
# print n_values
# print m
print 'Number of runs in file ', N_values

if 'slice_start' not in dir():
  slice_start = 0

if 'slice_end' not in dir() or slice_end == 0:
  slice_end = N_values


#Define which files need tweaking
code_files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
         'math', 'meshgen', 'micro']
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

if task == 'fit':
  hout = open('gnuplot_scripts/fitTaylorDF.run.plt', 'w')
# Fortran uses d to indicate double precision in scientific notation:
# eg, 3.14d0. This converts these floating point numbers to using an e
# instead of a d, when called by re.sub.
  def convert_d_to_e(matchobj):
    return matchobj.group().replace('d','e')

if task == 'timeseries':
  gnuplotcmd = "plot "

# This is the main loop, setting up each of the runs
for values in listruns.ItRunValues(list_values, tokens, n_values, N_values, m, 
                              pattern, tokendict, slice_start, slice_end):
# Do the string replace operations on the values themselves
  cd = values[tokendict['$label$']]
  wd = os.path.join('.', 'batch', cd)

  if task == 'run' or task == 'rerun':
    print cd
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
    if system == 'tara':
      hin = open('tara.serial.raw.slurm','r')
      houtcode = open(os.path.join(wd, 'tara.serial.run.slurm'), 'w')
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

    if system == 'popledebug':
      hin = open('qsub_pople_debug.raw','r')
      houtcode = open(os.path.join(wd, 'qsub_pople_debug.run'), 'w')
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
      fortran_command = 'gfortran -Wall -g -m64 -o cell' + cd
    elif system == 'gfortranopenmp':
      fortran_command = 'gfortran -Wall -fopenmp -g -o cell' + cd
    elif system == 'sun':
      fortran_command = 'f95 -o cell' + cd
    elif system == 'ifort':
      fortran_command = 'ifort -o cell' + cd
    elif system == 'hpc':
      fortran_command = 'mpif90 -o cell' + cd
    elif system == 'tara':
      fortran_command = 'gfortran -o cell' + cd
    elif system == 'pople' or system == 'popledebug':
      fortran_command = 'ifort -o cell' + cd
    for file in code_files:
        fortran_command = fortran_command + ' ' + file + file_out_suffix
    os.system(fortran_command)
    if system == 'gfortran' or system == 'gfortranopenmp' or system == 'ifort'\
           or system == 'sun':
      os.system('./cell' + cd)
    elif system == 'hpc':
      os.system('qsub qsub_script.run')
    elif system == 'tara':
      os.system('sbatch tara.serial.run.slurm')
    elif system == 'pople':
      os.system('qsub qsub_pople.run')
    elif system == 'popledebug':
      os.system('qsub qsub_pople_debug.run')
    os.chdir(os.path.join('..', '..'))

  elif task == 'extract':
    print os.path.join(wd, extractfile)
    if os.path.exists(os.path.join(wd, extractfile)):
      hin = open(os.path.join(wd, extractfile),'r')
      hout.write(wd + ' ' + hin.readlines()[-1])
      hin.close()
  elif task == 'plot':
    print os.path.join(wd, extractfile)
    if os.path.exists(os.path.join(wd, extractfile)):
      hin = open(os.path.join(wd, extractfile),'r')
      r = values[tokendict['$r$']]
      capillary_no = values[tokendict['$capillary_no$']]
      line = r + ' ' + capillary_no + hin.readlines()[-1]
      hout.write(wd + ' ' + line)
      hin.close()
  elif task == 'list':
    print cd
  elif task == 'foreach':
    print cd
    os.chdir(wd)
    os.system(options.foreach)
    os.chdir(os.path.join('..', '..'))
  elif task == 'fit':
    print cd
    hin = open('gnuplot_scripts/fitTaylorDF.raw.plt')
    for line in hin.readlines():
      for j in range(0,len(tokens)):
        line = line.replace(tokens[j], values[j])
      line = re.sub(r'\d.?d-?\d', convert_d_to_e, line)
# 
      line = line.replace('$cd$', cd)
      line = line.replace('$wd$', wd)
      hout.write(line)
    hin.close()
  elif task == 'clone':
    print cd
    if not os.path.exists(os.path.join('Alex', 'clone', cd)):
      os.mkdir(os.path.join('Alex', 'clone', cd))
    os.chdir(wd)
    os.system('cp ' + clonefile + ' ' +
              os.path.join('..', '..', 'Alex', 'clone', cd))
    os.chdir(os.path.join('..', '..'))
  elif task == 'backup':
    print cd
    if os.path.exists(os.path.join('Alex', 'backup', cd)):
      os.remove(os.path.join('Alex', 'backup', cd))
    os.system('cp -R ' + os.path.join('batch', cd) + ' ' +
              os.path.join('Alex', 'backup'))
  elif task == 'restore':
    print cd
    if not os.path.exists(os.path.join('batch', cd)):
      os.system('cp -R ' + os.path.join('Alex', 'backup', cd) + ' ' +
                os.path.join('batch', cd))
    else:
      print 'Error: batch directory ' + cd
      print 'already exists, and will not be overwritten by the backup.'
      print 'Manually remove ' + os.path.join('batch', cd)
      print 'and try again.'
  elif task == 'timeseries':
    gnuplotcmd += "'" + os.path.join(wd, 'TaylorDF__00001.txt') + "' w l, "

if task == 'extract' or task == 'plot':
  hout.close()
  hin = open('extract' + extractfile, 'r')
  print hin.read()
  hin.close()

if task == 'fit':
  hout.close()
  os.system('gnuplot gnuplot_scripts/fitTaylorDF.run.plt')

if task == 'timeseries':
#  hout = open('gnuplot_scripts/timeseries-temp-' + myfile + '.txt', 'w')
  hout = open('gnuplot_scripts/timeseries-temp.txt', 'w')
  hout.write(gnuplotcmd[:-2] + '\n')
  hout.write("pause -1 'Hit return to continue'")
  hout.close()
  os.system('gnuplot gnuplot_scripts/timeseries-temp.txt')
