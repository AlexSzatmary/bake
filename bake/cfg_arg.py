#!/usr/bin/env python
"""
This module handles importing options from files and arguments from the 
command line.
"""

import argparse
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


def make_argparser(argv):
    argparser = argparse.ArgumentParser(prog="bake", epilog=__doc__)
    # Core tasks
    argparser.add_argument('--file', '-f',
                         help="""Bake parameter file to operate from""")
    argparser.add_argument(
        '--bake_file', '-b', action='append',
        help="""File to bake. This can be a filename or a glob. This overrides
        the list bake_files in bake.cfg if it is present. This option can be
        used repeatedly, with each -b flag specifying another file or glob to
        bake.""")

    argparser.add_argument('--mix', '-m',
                         help="""Mix parameters into code files.""",
                         action='store_true')
    argparser.add_argument(
      '--remix', '-M',
      help="""Like mix, but does not make new directories.""",
      action='store_true')
    argparser.add_argument('--list', '-l', action='store_true',
                         help="""Lists the jobs that would be operated on with
                         the given parameter file and options.""")
    # Modifiers
    argparser.add_argument(
        '--overwrite', '-o', action='append',
        help="""Overwrite a line in a batch parameter file, eg,"-o
'\$foo\$;bar;baz'" replaces a parameter line starting with"$foo$" with
"$foo$;bar;baz". This option can be used repeatedly.  (Note: if the parameter
specified is absent from the file, the new line will simply be added to the
options in the file, it won't overwrite anything.)
        """)
    argparser.add_argument(
        '--slice', '-s', help="""Selects a subset of the runs specified in the
file, eg, -s 5-9 does runs 5, 6, 7, and 8 out of however many runs would be
referred to in the given file.
        """)
    argparser.add_argument('--execute', '-e', help="""Execute a command in each
parameter specified, eg, "tail foo.txt"
        """)
#    argparser.add_argument('--bake_files', help="""Files to mix""")
    return argparser


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


