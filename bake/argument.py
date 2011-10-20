#!/usr/bin/env python
"""
This module handles importing options from files and arguments from the 
command line.
"""

import argparse
import bakedefaults


class BakeArgparser(argparse.ArgumentParser):
    def __init__(self):
        argparse.ArgumentParser.__init__(self, prog="bake", epilog=__doc__,
                                        fromfile_prefix_chars='+')
        # Core tasks
        self.add_argument('--file', '-f',
                             help="""Bake parameter file to operate from""")
        self.add_argument(
            '--bake_file', '-b', action='append',
            help="""File or comma-separated list of files to bake. This can be
            a filename or a glob. This option can be used repeatedly, with
            each -b flag specifying another file or glob to bake.""")
        self.add_argument('--prefix', '-p', help="""The directory with the grid
            of jobs.""", default='batch')
        self.add_argument('--mix', '-m',
                             help="""Mix parameters into code files.""",
                             action='store_true')
        self.add_argument(
          '--remix', '-M',
          help="""Like mix, but does not make new directories.""",
          action='store_true')
        self.add_argument('--list', '-l', action='store_true',
                             help="""Lists the jobs that would be operated on
                             with the given options.""")
        # Modifiers
        self.add_argument(
            '--overwrite', '-o', action='append',
            help="""Overwrite a line in a batch parameter file, eg,"-o
            '@foo@;bar;baz'" replaces a parameter line starting with"@foo@"
            with "$foo$;bar;baz". This option can be used repeatedly.  (Note:
            if the parameter specified is absent from the file, the new line
            will simply be added to the options in the file, it won't overwrite
            anything.)
            """)
        self.add_argument(
            '--slice', '-s', help="""Selects a subset of the runs specified in
            the file, eg, -s 5-9 does runs 5, 6, 7, and 8 out of however many
            runs would be referred to in the given file.
            """)
        self.add_argument('--execute', '-e', help="""Execute a command in each
            parameter specified, eg, "tail foo.txt"
            """)
        self.add_argument('--key_start', default='@')
        self.add_argument('--key_end', default='@')

    def convert_arg_line_to_args(self, arg_line):
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg


def process_arguments(args):
    """
    This is used to turn the args object into something that can be used
    better later in the code. It does error checking and massaging of the
    command line arguments. It checks to make sure a file was input, for
    example.
    """
    try:
        # A file must be selected to operate on
#         if not args.file:
#             raise Exception('No batch parameter file specified')
        # Perform operation on a Slice of the runs
        if args.slice:
            if '-' in args.slice:
                (args.slice_start, args.slice_end) = \
                    args.slice.split('-')
            else:
                args.slice_start = 0
                args.slice_end = int(args.slice)
            if args.slice_start == '':
                args.slice_start = 0
            if args.slice_end == '':
                args.slice_end = 0
            args.slice_start = int(args.slice_start)
            args.slice_end = int(args.slice_end)
        else:
            args.slice_start = 0
            args.slice_end = 0

        if args.remix:
            if args.mix:
                print('You set both the mix and remix flags; to remix, you '
                      'just need to set the remix flag.')
            args.mix = True
        
        if args.bake_file:
            bake_file = []
            for a in args.bake_file:
                bake_file.extend(a.split(','))
            args.bake_file = bake_file

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


