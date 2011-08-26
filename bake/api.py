#!/usr/bin/python
# api.py

"""
This module takes bake parameter files, which are lists of tags and values,
then does something for each combination of values.
Okay, that's pretty vague and broad.
It's useful, honest, for doing repetitive find-and-replace operations,
wrangling data out of oodles of really-similar-but-sublty-different-in-two-
variables-but-not-the-other-five sets of data, doing the accounting on
submitting jobs to all of the different kinds of supercomputers in your life,
and making plots of y vs x for a given z, but then b vs t for a given y.

It's like a little robot that does repetitive things for you so you don't get
tenosynovitis.
"""

import os
import os.path
import mix
import optparse
import re
import glob
import load
import bakedefaults

def load_config():
    """
    Reads bake.cfg and makes a dict out of it
    """
    import ConfigParser
    c = ConfigParser.SafeConfigParser()
    c.read(bakedefaults.CFGFILE)
    d = {}
    for section in c.sections():
        d[section] = dict(c.items(section))
        for k in d[section]:
            # Handle newlines in values in a more useful way
            # Make comma-separated lists into list objects
            d[section][k] = d[section][k].replace('\n', '').split(',')
    return d


def make_optparser():
    optparser = optparse.OptionParser(prog="bake", epilog = __doc__)
    # Core tasks
    optparser.add_option('--file', '-f',
                         help="""Bake parameter file to operate from""")
    optparser.add_option(
        '--bake_file', '-b', action='append',
        help="""File to bake. This can be a filename or a glob. This overrides
        the list bake_files in bake.cfg if it is present. This option can be
        used repeatedly, with each -b flag specifying another file or glob to
        bake..""")

    optparser.add_option('--mix', '-m',
                         help="""Mix parameters into code files.""",
                         action='store_true')
    optparser.add_option(
      '--remix', '-M',
      help="""Like mix, but does not make new directories.""",
      action='store_true')
    optparser.add_option('--list', '-l', action='store_true',
                         help="""Lists the jobs that would be operated on with
                         the given parameter file and options.""")
    # Modifiers
    optparser.add_option(
        '--overwrite', '-o', action='append',
        help="""Overwrite a line in a batch parameter file, eg,"-o
'\$foo\$;bar;baz'" replaces a parameter line starting with"$foo$" with
"$foo$;bar;baz". This option can be used repeatedly.  (Note: if the parameter
specified is absent from the file, the new line will simply be added to the
options in the file, it won't overwrite anything.)
        """)
    optparser.add_option(
        '--slice', '-s', help="""Selects a subset of the runs specified in the
file, eg, -s 5-9 does runs 5, 6, 7, and 8 out of however many runs would be
referred to in the given file.
        """)
    optparser.add_option('--execute', '-e', help="""Execute a command in each
parameter specified, eg, "tail foo.txt"
        """)
    return optparser


def process_options(options):
    """
    This is used to turn the options object into something that can be used
    better later in the code. It does error checking and massaging of the
    command line arguments. It checks to make sure a file was input, for
    example.
    """
    try:
        # A file must be selected to operate on
#         if not options.file:
#             raise Exception('No batch parameter file specified')
        # Perform operation on a Slice of the runs
        if options.slice:
            if '-' in options.slice:
                (options.slice_start, options.slice_end) = \
                    options.slice.split('-')
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

        if options.remix:
            if options.mix:
                print('You set both the mix and remix flags; to remix, you '
                      'just need to set the remix flag.')
            options.mix = True

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


def make_grid(config, options, lines):
    """
    This is the interface between the internals of bake.mix and what people
    would practically use.
    """
    grid = mix.Grid(lines)
    grid.set_slice(options.slice_start, options.slice_end)
    grid.set_key_pattern(config, options.file)
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
