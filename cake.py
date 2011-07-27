#!/usr/bin/python
# cake.py
# Alex Szatmary
# This is bake for the cell code: 'cell'[0] + 'bake'[1:] = cake

import sys
sys.path.insert(0, '.')
sys.path.insert(0, './Alex')
import trim_tail
import os
import os.path
import re
import optparse
import mix
import projectPrefs
import bake
print os.getcwd()

def cake():
  optparser = bake.make_optparser()
  projectPrefs.set_opt_parser(optparser)
  optparser.add_option('--run', '-r',
                       help="""Start a run; specify a system to run on.
                       eg, "-r foo" asks to run the code on system foo.
                       A system must be specified. Currently used systems
                       include "hpc", "pople", "tara", and "sun".""")
  optparser.add_option('--rerun', '-R',
                       help="""Like -r, but restarts a stopped run.
                       The only difference between this and -r is that new
                       directories are not made.""",)
  optparser.add_option('--trimtail', '-T', action='store_true',
                       help="Deletes all datapoints after the most recent "
                       "checkpoint.")

  options, arguments = optparser.parse_args()

  bake.process_options(options)

  ## Prefs
  config = projectPrefs.InitializeOptions
  bake.load_config(config)
  task = ''

  try:
    if options.run or options.rerun:
      if task:
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

    if options.trimtail:
      if task:
        raise Exception('Multiple tasks requested')
      else:
        task = 'trimtail'

  except Exception, data:
    if data[0] == 'Invalid system specified':
      print data[0]
      exit(-1)
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

  print(options)
  if options.overwrite:
    lines = options.overwrite
  else:
    lines = []

  # in lines
  hin = open(options.file,'r')
  lines += hin.readlines()
  hin.close()

  (label, tokens, 
   mixIterator) = bake.make_iterator(config.label_tag, 
                                     re.compile(config.pattern), 
                                     lines, options.slice_start, 
                                     options.slice_end)
   

  ## This is the main loop, iterating over each set of values
  if options.run or options.rerun:
    for values in mixIterator:
    # Do the string replace operations on the values themselves
      cd = values[label]
      wd = os.path.join('.', 'batch', cd)

      if task == 'run' or task == 'rerun' or task == 'mix':
        if task == 'run' or task == 'mix':
          os.mkdir(wd)
      # String replace the tokens for the values
        for f in config.code_files:
            hin = open(f + config.file_in_suffix,'r')
            houtcode = open(os.path.join(wd, f + config.file_out_suffix), 'w')
            for line in hin.readlines():
                for j in range(0,len(tokens)):
                    line = line.replace(tokens[j], values[j])
                houtcode.write(line)
            hin.close()
            houtcode.close()
        for f in config.visual_files:
            hin = open(f + config.visual_file_in_suffix,'r')
            houtcode = open(os.path.join(wd, f + config.visual_file_out_suffix), 'w')
            for line in hin.readlines():
                for j in range(0,len(tokens)):
                    line = line.replace(tokens[j], values[j])
                houtcode.write(line)
            hin.close()
            houtcode.close()
        if task == 'run' or task == 'rerun':
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
          for f in config.code_files:
              fortran_command = (fortran_command + ' ' + f + 
                                 config.file_out_suffix)
          os.system(fortran_command)
          if system == 'gfortran' or system == 'gfortranopenmp' or \
                system == 'ifort' or system == 'sun':
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
  elif task == 'trimtail':
    for values in mixIterator:
    # Do the string replace operations on the values themselves
      cd = values[label]
      wd = os.path.join('.', 'batch', cd)
      trim_tail.trim_tail_directory(wd)
  else:
    bake.default_loop(label, tokens, mixIterator, config, options)
  
if __name__ == '__main__':
  cake()
