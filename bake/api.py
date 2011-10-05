#!/usr/bin/python
# api.py

"""
This module takes bake parameter files, which are lists of tags and values,
then does something for each combination of values.
Okay, that's pretty vague and broad.
It's useful, honest, for doing repetitive find-and-replace operations,
wrangling data out of oodles of really-similar-but-subtly-different-in-two-
variables-but-not-the-other-five sets of data, doing the accounting on
submitting jobs to all of the different kinds of supercomputers in your life,
and making plots of y vs x for a given z, but then b vs t for a given y.

It's like a little robot that does repetitive things for you so you don't get
tenosynovitis.
"""

import os
import os.path
import mix
import re
import glob
import load
import cfg_arg

def make_grid(lines, options, config):
    """
    This is the interface between the internals of bake.mix and what people
    would practically use.
    """
    grid = mix.Grid(lines)
    if options:
        grid.set_slice(options.slice_start, options.slice_end)
    else:
        grid.set_slice()
    if options:
        grid.set_key_pattern(config, options.file)
    else:
        grid.set_key_pattern(config)
    return grid


def default_make_grid(lines):
    """
    New goal: make this function useless
    """
    grid = mix.Grid(lines)
    grid.set_slice()
    grid.set_key_pattern()
    return grid


def default_loop(grid, options):
    """
    This does a loop over each item in mixIterator, that is, each combination
    of the possible values.
    """
    for values in grid.mix_iterator():
    # Do the string replace operations on the values themselves
        cd = grid.get_label()
        wd = os.path.join('.', 'batch', cd)

        if options.list:
            print(cd)
        # mix is the main special thing done by this loop:
        # it takes the files, does a wild find and replace on them to take out
        # the keys and swap in corresponding values, and spits the files out
        # in ./batch/ somewhere.
        if options.mix:
            if not options.remix:
                os.mkdir(wd)
            # String replace the keys for the values
            if options.bake_file:
                for g in options.bake_file:
                    for f in glob.glob(g):
                        hin = open(f + options.file_in_suffix, 'r')
                        houtcode = open(
                            os.path.join(wd, f + options.file_out_suffix), 'w')
                        for line in hin.readlines():
                            line = grid.replace(line)
                            houtcode.write(line)
                        hin.close()
                        houtcode.close()
        if options.execute:
            os.chdir(wd)
            os.system(grid.replace(options.execute))
            os.chdir(os.path.join('..', '..'))
