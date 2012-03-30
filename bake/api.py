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
import argument

def make_grid(lines, args):
    """
    This is the interface between the internals of bake.mix and what people
    would practically use.
    """
    grid = mix.Grid(lines)
    if args:
        grid.set_slice(args.slice_start, args.slice_end)
    else:
        grid.set_slice()
    grid.set_key_pattern(args.key_start, args.key_end, args.file)
    return grid


def default_make_grid(lines):
    """
    New goal: make this function useless
    """
    grid = mix.Grid(lines)
    grid.set_slice()
    grid.set_key_pattern()
    return grid


def default_loop(grid, args):
    """
    This does a loop over each item in mixIterator, that is, each combination
    of the possible values.
    """
    # Optionally record a manifest about this job, what the effective bp file
    # was
    if args.grid_manifest:
        cd = os.getcwd()
        os.chdir(args.grid_manifest_dir)
        grid.write_manifest()
        os.chdir(cd)

    for values in grid.mix_iterator():
    # Do the string replace operations on the values themselves
        cd = grid.get_label()
        if not os.path.exists(args.prefix):
            raise OSError, 'Bake requires the directory, ' + args.prefix
        wd = os.path.join(args.prefix, cd)

        if args.list:
            print(cd)
        # mix is the main special thing done by this loop:
        # it takes the files, does a wild find and replace on them to take out
        # the keys and swap in corresponding values, and spits the files out
        # in ./batch/ somewhere.
        if args.mix:
            if not args.remix:
                os.mkdir(wd)
            # String replace the keys for the values
            if args.bake_file:
                for g in args.bake_file:
                    for f in glob.glob(g):
                        hin = open(f, 'r')
                        houtcode = open(os.path.join(wd, f), 'w')
                        for line in hin.readlines():
                            line = grid.replace(line)
                            houtcode.write(line)
                        hin.close()
                        houtcode.close()
        if args.execute:
            os.chdir(wd)
            for e in args.execute:
                os.system(grid.replace(e))
            os.chdir(os.path.join('..', '..'))
        # Optionally record a bp file representing a single job for each job
        if args.job_manifest:
            os.chdir(wd)
            grid.write_job_manifest()
            os.chdir(os.path.join('..', '..'))
            
