#!/usr/bin/python
# cmdline.py
# This is a custom front-end for bake for the poisson example.

import sys
import os
import os.path
import argparse
sys.path.insert(0, '../../../bake')
import bake
import bake.mix


def main(argv=sys.argv[1:]):
    """
    Main routine for cmdline.py
    """

    # Generate an optparser object with the bake default options
    argparser = bake.argument.BakeArgparser()
    if os.path.exists('bake.arg'):
        argv = ['+bake.arg'] + argv

    # Add my own options
    argparser.add_argument(
        '--poisson', '-P', action='store_true', help="""
        Run the poisson solver
        """
        )
    argparser.add_argument(
        '--compare_ideal', '-c', action='store_true', help="""
        Compare the numerical solution with the ideal value, giving the
        L_infinity norm for each job specified in the bp file; this data is
        exported as a table in the error.txt file.
        """
        )

    args = argparser.parse_args(argv)

    # Bake tweaks the options object, it checks for errors, does
    # arithmetic for slice, and so on
    bake.argument.process_arguments(args)

    task = ''

    # Exception handling is mostly to make sure that multiple, conflicting,
    # tasks aren't requested. For example, a regular execute plus the poisson
    # task can't both be requested.
    try:
        if args.poisson:
            if task:
                raise Exception('Multiple tasks requested')
            # This makes the poisson option act like a simple execute option,
            # but with a built-in command
            task = 'execute'
            args.execute = ['python poisson.py']
            # Uncomment the below lines to force bake to always mix or list
            # args.mix = True
            # args.list = True
        if args.compare_ideal:
            if task:
                raise Exception('Multiple tasks requested')
            error_table = []
            hout = open('error.txt', 'w')

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
    ## End processing of command line parameters

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

    # If you were to have many custom commands, set this up like a series of
    # elif's
    if args.compare_ideal:
        # To loop over the grid of jobs, do a for loop on mixIterator
        for values in grid.mix_iterator():
            cd = grid.get_label()
            wd = os.path.join('.', 'batch', cd)

            if args.list:
                print(cd)
            os.chdir(wd)
            h = os.popen('python compare_ideal.py')
            error = h.read()
            h.close()
            os.chdir('../..')

            table_line = '@n@ ' + error
            # This for loop is set up so that any tag could be put in
            # the default table_line and replaced for its value
            table_line = grid.replace(table_line)
            error_table.append(table_line)
    # You can add other tasks like this
    # elif args.other_task:
    #     code for other task
    # bake.default_loop preserves bake's regular behavior plus this is how
    # the poisson task is run
    else:
        bake.default_loop(grid, args)

    # Post-processing
    if args.compare_ideal:
        hout = open('error.txt', 'w')
        hout.writelines(error_table)
        hout.close()

if __name__ == '__main__':
    main()
