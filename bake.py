#!/usr/bin/python
# bake.py
# Alex Szatmary
# This script takes bake parameter files, which are lists of tags and values,
# then does something for each combination of values.
# Okay, that's pretty vague and broad.
# It's useful, honest, for doing repetitive find-and-replace operations,
# wrangling data out of oodles of really-similar-but-sublty-different-in-two-
# variables-but-not-the-other-five sets of data, doing the accounting on
# submitting jobs to all of the different kinds of supercomputers in your life,
# and making plots of y vs x for a given z, but then b vs t for a given y.
#
# It's like a little robot that does repetitive things for you so you don't get
# carpal tunnel.

import sys
sys.path.insert(0, '.')
import os
import os.path
import re
import mix
import optparse
import projectPrefs

def bake():
  optparser = optparse.OptionParser()
  # Core tasks 
  optparser.add_option('--mix', '-m',
                       help="""Mix parameters into code files.""", 
                       action='store_true')
  optparser.add_option('--list', '-l', action='store_true',
                       help="""Lists the jobs that would be operated on with the
                       given parameter file and options.""")
  # Modifiers
  optparser.add_option('--overwrite', '-o', action='append',
                       help="""Overwrite a line in a batch parameter file,
                       eg, "-o '\$foo\$;bar;baz'" replaces a parameter line
                       starting with "$foo$" with "$foo$;bar;baz".
                       This option can be used repeatedly.
                       (Note: if the parameter specified is absent from the file,
                       the new line will simply be added to the options in the
                       file, it won't overwrite anything.)
                       """)
  optparser.add_option('--slice', '-s',
                       help="""Selects a subset of the runs specified in the
                       file, eg, -s 5:9 does runs 5, 6, 7, and 8 out of however
                       many runs would be referred to in the given file.""")

  # Cultural tasks
  optparser.add_option('--foreach', '-E',
                       help="""Execute a command in each job specified, eg, "tail
                       TaylorDF__0001.txt"
                       """)
  optparser.add_option('--backup', '-b', action='store_true',
                       help="""Backs up selected runs to Alex/backup. This is
                       helpful to do between runs, after checking to make sure
                       that the checkpoint and data files are saved safely.
                       """)
  optparser.add_option('--restore', '-t', action='store_true',
                       help="""Brings back runs backed up by the backup option.
                       """)
  projectPrefs.set_opt_parser(optparser)
  options, arguments = optparser.parse_args()

  #todo Separate out the core from the prefs
  try:
    if options.mix:
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      if options.mix:
        task = 'mix'

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

    if options.foreach:
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      task = 'foreach'

    if options.backup:
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      task = 'backup'

    if options.restore:
      if 'task' in dir():
        raise Exception('Multiple tasks requested')
      task = 'restore'

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


  ## End processing of command line parameters
  ## Prepare for big loop

  if options.overwrite:
    lines = options.overwrite
  else:
    lines = []

  hin = open(myfile,'r')
  lines += hin.readlines()
  hin.close()

  (tokens, list_values, n_values, N_values, tokendict, m) = \
      mix.parseBPlines(lines)

  if 'slice_start' not in dir():
    slice_start = 0

  if 'slice_end' not in dir() or slice_end == 0:
    slice_end = N_values


  ## Prefs
  if task != 'foreach':
    print 'Number of runs in file ', N_values

  opts = projectPrefs.InitializeOptions


  ## This is the main loop, iterating over each set of values
  #todo Separate out the core from the prefs
  for values in mix.ItRunValues(list_values, tokens, n_values, N_values, m, 
                                opts.pattern, tokendict, slice_start, slice_end):
  # Do the string replace operations on the values themselves
    cd = values[tokendict[opts.label_tag]]
    wd = os.path.join('.', 'batch', cd)

    if task == 'run' or task == 'rerun' or task == 'mix':
      print cd
      if task == 'run' or task == 'mix':
        os.mkdir(wd)
    # String replace the tokens for the values
      for f in opts.code_files:
          hin = open(f + opts.file_in_suffix,'r')
          houtcode = open(os.path.join(wd, f + opts.file_out_suffix), 'w')
          for line in hin.readlines():
              for j in range(0,len(tokens)):
                  line = line.replace(tokens[j], values[j])
              houtcode.write(line)
          hin.close()
          houtcode.close()
    elif task == 'list':
      print cd
    elif task == 'foreach':
      #print cd
      os.chdir(wd)
      os.system(options.foreach)
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

  ## Post-loop
  #todo Separate out the core from the prefs

  if task == 'extract':
    hout.close()
    hin = open('extract' + extractfile, 'r')
    print hin.read()
    hin.close()

if __name__ == '__main__':
  bake()
