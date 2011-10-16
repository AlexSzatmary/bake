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
import os.path

def main(argv=sys.argv[1:]):
    # Set up command line argument options
    argparser = bake.argument.BakeArgparser()
    if os.path.exists('bake.arg'):
        argv = ['+bake.arg'] + argv
    args = argparser.parse_args(argv)
    bake.argument.process_arguments(args)

    ## End processing of command line parameters
    ## Prepare for big loop

    # Load bake parameter file
    if args.file:
        lines = bake.load.load_file(args.file)
    else:
        lines = []

    # The overwrite command pushes lines onto the top of the bake parameter
    # file
    if args.overwrite:
        lines.extend(bake.load.load(l for l in args.overwrite))

    # This grid object is kind of the core of bake.
    grid = bake.make_grid(lines, args)

    ## This is the main loop, iterating over each set of values
    bake.default_loop(grid, args)

if __name__ == '__main__':
    main()
