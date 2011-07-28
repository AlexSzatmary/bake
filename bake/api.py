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

import os
import os.path
import mix
import optparse
import re

class Config:
  pass

def load_config():
  import ConfigParser
  c = ConfigParser.SafeConfigParser()
  c.read('bake.cfg')
  config = Config()
  config.label_tag = c.get('label', 'label_tag')
  config.pattern = c.get('label', 'pattern')

  config.code_files = c.get('filenames', 'code_files').split(',')
#  print(config.code_files)
  config.file_in_suffix = c.get('filenames', 'file_in_suffix')
#  print(config.file_in_suffix)
  config.file_out_suffix = c.get('filenames', 'file_out_suffix')
  return config

def make_optparser():
  optparser = optparse.OptionParser()
  # Core tasks 
  optparser.add_option('--file', '-f',
                       help="""Bake parameter file to operate from""")
  optparser.add_option('--mix', '-m',
                       help="""Mix parameters into code files.""", 
                       action='store_true')
  optparser.add_option('--list', '-l', action='store_true',
                       help="""Lists the jobs that would be operated on with
                       the given parameter file and options.""")
  # Modifiers
  optparser.add_option('--overwrite', '-o', action='append',
                       help="""Overwrite a line in a batch parameter file, eg,
                       "-o '\$foo\$;bar;baz'" replaces a parameter line
                       starting with "$foo$" with "$foo$;bar;baz".  This option
                       can be used repeatedly.  (Note: if the parameter
                       specified is absent from the file, the new line will
                       simply be added to the options in the file, it won't
                       overwrite anything.)
                       """)
  optparser.add_option('--slice', '-s',
                       help="""Selects a subset of the runs specified in the
                       file, eg, -s 5-9 does runs 5, 6, 7, and 8 out of however
                       many runs would be referred to in the given file.""")
  optparser.add_option('--execute', '-e',
                       help="""Execute a command in each job specified, eg,
                       "tail TaylorDF__0001.txt"
                       """)

#   # Cultural tasks
#   optparser.add_option('--backup', '-b', action='store_true',
#                        help="""Backs up selected runs to Alex/backup. This is
#                        helpful to do between runs, after checking to make sure
#                        that the checkpoint and data files are saved safely.
#                        """)
#   optparser.add_option('--restore', '-t', action='store_true',
#                        help="""Brings back runs backed up by the backup option.
#                        """)
  return optparser

def process_options(options):
  try:
    if not options.file:
      raise Exception('No batch parameter file specified')
    # Perform operation on a Slice of the runs      
    if options.slice:
      if '-' in options.slice:
        (options.slice_start, options.slice_end) = options.slice.split('-')
      else:
        options.slice_start = 0
        options.slice_end = int(options.slice)
      if options.slice_start == '':
        options.slice_start = 0
      if options.slice_end == '':
        options.slice_end = 0
      options.slice_start = int(options.slice_start)
      options.slice_end = int(options.slice_end)
    else:
      options.slice_start = 0
      options.slice_end = 0

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

def make_iterator(label_tag, pattern, lines, slice_start, slice_end):
  (tokens, list_values, n_values, N_values, tokendict) = \
      mix.parseBPlines(lines)
  mixIterator = mix.ItRunValues(list_values, tokens, n_values, N_values, 
                                pattern, tokendict, slice_start, slice_end)
  return (tokendict[label_tag], tokens, mixIterator)

def default_loop(label, tokens, mixIterator, config, options):
  for values in mixIterator:
  # Do the string replace operations on the values themselves
    cd = values[label]
    wd = os.path.join('.', 'batch', cd)

    if options.list:
      print(cd)
    if options.mix:
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
    if options.execute:
      os.chdir(wd)
      os.system(options.execute)
      os.chdir(os.path.join('..', '..'))

# def bake():
#   optparser = make_optparser()
#   options, arguments = optparser.parse_args()

#   process_options(options)

#   ## Prefs
#   load_config(config)

#   ## End processing of command line parameters
#   ## Prepare for big loop

#   if options.overwrite:
#     lines = options.overwrite
#   else:
#     lines = []

#   # in lines
#   hin = open(options.file,'r')
#   lines += hin.readlines()
#   hin.close()

#   (label, tokens, 
#    mixIterator) = make_iterator(config.label_tag, re.compile(config.pattern), 
#                                      lines, options.slice_start, 
#                                      options.slice_end)
   
#   ## This is the main loop, iterating over each set of values
#   default_loop(label, tokens, mixIterator, config, options)

#     if options.backup:
#       if os.path.exists(os.path.join('Alex', 'backup', cd)):
#         os.remove(os.path.join('Alex', 'backup', cd))
#       os.system('cp -R ' + os.path.join('batch', cd) + ' ' +
#                 os.path.join('Alex', 'backup'))
#     elif options.restore:
#       if not os.path.exists(os.path.join('batch', cd)):
#         os.system('cp -R ' + os.path.join('Alex', 'backup', cd) + ' ' +
#                   os.path.join('batch', cd))
#       else:
#         print 'Error: batch directory ' + cd
#         print 'already exists, and will not be overwritten by the backup.'
#         print 'Manually remove ' + os.path.join('batch', cd)
#         print 'and try again.'

#if __name__ == '__main__':
#  bake()
