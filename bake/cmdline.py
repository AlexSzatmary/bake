#!/usr/bin/env python
# encoding: utf-8
"""
This is the command line interface for bake. For people who want to take
bake.py and extend it for their own circumstances, modifying the main routine
in this module is probably the best place to start.
"""

import api as bake
import load
import sys


def main(args=sys.argv[1:]):
    # Set up command line argument options
    optparser = bake.make_optparser()
    options, arguments = optparser.parse_args()

    bake.process_options(options)

    ## Configuration is stored in the bake.cfg file in the current directory
    config = bake.load_config()

    ## End processing of command line parameters
    ## Prepare for big loop

    # The overwrite command pushes lines onto the top of the bake parameter
    # file
    if options.overwrite:
        lines = options.overwrite
    else:
        lines = []

    # Load bake parameter file
    load.load(l for l in lines)
    lines.extend(load.load_file(options.file))

    # This mixIterator object is kind of the core of bake.
    (label, tokens,
     mixIterator) = bake.make_iterator(config['label']['label_tag'],
                                       config['label']['pattern'],
                                       lines, options.slice_start,
                                       options.slice_end)

    ## This is the main loop, iterating over each set of values
    bake.default_loop(label, tokens, mixIterator, config, options)

if __name__ == '__main__':
    main()
