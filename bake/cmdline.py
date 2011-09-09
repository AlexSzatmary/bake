#!/usr/bin/env python
# encoding: utf-8
"""
This is the command line interface for bake. For people who want to take
bake.py and extend it for their own circumstances, modifying the main routine
in this module is probably the best place to start.
"""

import argparse
# If you make your own bake frontend, just do
# import bake
import api as bake
# and do bake.load
# import load
import sys

def main(argv=sys.argv[1:]):
    # Set up command line argument options
#    options = bake.cfg_arg.get_options(['@bake.cfg'] + argv)
    argparser = bake.cfg_arg.make_argparser(argv)
    args = argparser.parse_args(argv)
    bake.cfg_arg.process_options(args)

    ## Configuration is stored in the bake.cfg file in the current directory
    config = bake.cfg_arg.load_config()

    ## End processing of command line parameters
    ## Prepare for big loop

    # Load bake parameter file
    if args.file:
        lines = bake.load.load_file(args.file)
    else:
        lines = []

    if not args.bake_file:
        if 'filenames' in config and 'bake_files' in config['filenames']:
            args.bake_file = config['filenames']['bake_files']

    #warn This adds secret options to the options object.
    args.file_in_suffix = ''
    args.file_out_suffix = ''
    if 'filenames' in config and 'file_in_suffix' in config['filenames']:
        args.file_in_suffix = config['filenames']['file_in_suffix'][0]
    if 'filenames' in config and 'file_out_suffix' in config['filenames']:
        args.file_in_suffix = config['filenames']['file_out_suffix'][0]

    # The overwrite command pushes lines onto the top of the bake parameter
    # file
    if args.overwrite:
        lines.extend(bake.load.load(l for l in args.overwrite))

    # This mixIterator object is kind of the core of bake.
    grid = bake.make_grid(lines, args, config)

    ## This is the main loop, iterating over each set of values
    bake.default_loop(grid, args)

if __name__ == '__main__':
    main()
