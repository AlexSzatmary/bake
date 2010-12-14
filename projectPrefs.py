#!/usr/bin/python

def setoptparser(optparser):
    optparser.add_option('--extract', '-e',
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
