#!/usr/bin/env python
# encoding: utf-8

import api as bake
import sys

def main(args=sys.argv[1:]):
  optparser = bake.make_optparser()
  options, arguments = optparser.parse_args()

  bake.process_options(options)

  ## Prefs
  config = bake.load_config()

  ## End processing of command line parameters
  ## Prepare for big loop

  if options.overwrite:
    lines = options.overwrite
  else:
    lines = []

  # in lines
  hin = open(options.file,'r')
  lines += hin.readlines()
  hin.close()

  (label, tokens, 
   mixIterator) = bake.make_iterator(config['label']['label_tag'],
                                     config['label']['pattern'], 
                                     lines, options.slice_start, 
                                     options.slice_end)
   
  ## This is the main loop, iterating over each set of values
  bake.default_loop(label, tokens, mixIterator, config, options)

if __name__ == '__main__':
  main()
