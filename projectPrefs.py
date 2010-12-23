#!/usr/bin/python

import re, ConfigParser

def set_opt_parser(optparser):
    optparser.add_option('--run', '-r',
                         help="""Start a run; specify a system to run on.
                         eg, "-r foo" asks to run the code on system foo.
                         A system must be specified. Currently used systems
                         include "hpc", "pople", "tara", and "sun".""")
    optparser.add_option('--rerun', '-R',
                         help="""Like -r, but restarts a stopped run.
                         The only difference between this and -r is that new
                         directories are not made.""",)
    optparser.add_option('--extract', '-E',
			 help="""Extracts last line of the data file specified
			 here, for each run, prints these lines, and saves
                         them.
			 """)
    optparser.add_option('--plot', '-p',
			 help="Makes a datafile for a DF-m plot. "
			 "This is only useful right now for Alex's research, "
                         "but if you want code for something like this, let "
			 "him know.")
    optparser.add_option('--fit', '-f', action='store_true',
			 help="""Does a curve fit to 1-exp(t) using gnuplot for
			 Taylor DF.""")
    optparser.add_option('--timeseries', '-i', action='store_true',
			 help="""Make time-series DF plots for selected runs
			 """)

class InitializeOptions:
    # Define which files need tweaking
    code_files = ['cell', 'fluid', 'force', 'memb', 'rewr', 'visual', 'fvs',
                  'math', 'meshgen', 'micro']
    file_in_suffix = '.raw.f90'
    file_out_suffix = '.run.f90'
    label_tag = '$label$'
    pattern = re.compile('\$.+?\$')

    # These files are intended to make visualization with Matlab easier.
    visual_files = ['profilemovie']
    visual_file_in_suffix = '_raw.m'
    visual_file_out_suffix = '_run.m'
