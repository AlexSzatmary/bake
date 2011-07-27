#!/usr/bin/env python
# encoding: utf-8

import bake

def cmdline():
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
   mixIterator) = bake.make_iterator(config.label_tag, config.pattern, 
                                     lines, options.slice_start, 
                                     options.slice_end)
   
  ## This is the main loop, iterating over each set of values
  bake.default_loop(label, tokens, mixIterator, config, options)

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

if __name__ == '__main__':
  cmdline()
