#!/usr/bin/python
# cmdline.py
# This is a custom front-end for bake for the poisson example.

import sys
import os
import os.path
import optparse
import bake
import bake.mix


def main(args=sys.argv[1:]):
    """
    Main routine for cmdline.py
    """

    # Generate an optparser object with the bake default options
    optparser = bake.make_optparser()

    # Add my own options
    optparser.add_option(
        '--poisson', '-p', action='store_true', help="""
        Run the poisson solver
        """
        )
    optparser.add_option(
        '--compare_ideal', '-c', action='store_true', help="""
        Compare the numerical solution with the ideal value, giving the
        L_infinity norm for each job specified in the bp file; this data is
        exported as a table in the error.txt file.
        """
        )

    options, arguments = optparser.parse_args()

    # Bake tweaks the options object, it checks for errors, does
    # arithmetic for slice, and so on
    bake.process_options(options)

    # This loads the bake.cfg file as a dict, but processes some of the items.
    # For example, the bake_files string is converted to a list object
    config = bake.load_config()

    task = ''

    # Exception handling is mostly to make sure that multiple, conflicting,
    # tasks aren't requested. For example, a regular execute plus the poisson
    # task can't both be requested.
    try:
        if options.poisson:
            if task:
                raise Exception('Multiple tasks requested')
            # This makes the poisson option act like a simple execute option,
            # but with a built-in command
            task = 'execute'
            options.execute = ('python poisson.py')
            # Uncomment the below lines to force bake to always mix or list
            # options.mix = True
            # options.list = True
        if options.compare_ideal:
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
    if options.file:
        lines = bake.load.load_file(options.file)
    else:
        lines = []

    if not options.bake_file:
        if 'filenames' in config and 'bake_files' in config['filenames']:
            options.bake_file = config['filenames']['bake_files']

    #warn This adds secret options to the options object.
    options.file_in_suffix = ''
    options.file_out_suffix = ''
    if 'filenames' in config and 'file_in_suffix' in config['filenames']:
        options.file_in_suffix = config['filenames']['file_in_suffix'][0]
    if 'filenames' in config and 'file_out_suffix' in config['filenames']:
        options.file_in_suffix = config['filenames']['file_out_suffix'][0]

    # The overwrite command pushes lines onto the top of the bake parameter
    # file
    if options.overwrite:
        lines.extend(bake.load.load(l for l in options.overwrite))

    # This mixIterator object is kind of the core of bake.
    grid = bake.make_grid(lines, options, config)

    # If you were to have many custom commands, set this up like a series of
    # elif's
    if options.compare_ideal:
        # To loop over the grid of jobs, do a for loop on mixIterator
        for values in grid.mix_iterator():
            cd = grid.get_label()
            wd = os.path.join('.', 'batch', cd)

            if options.list:
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
    # elif options.other_task:
    #     code for other task
    # bake.default_loop preserves bake's regular behavior plus this is how
    # the poisson task is run
    else:
        bake.default_loop(grid, options)

    # Post-processing
    if options.compare_ideal:
        hout = open('error.txt', 'w')
        hout.writelines(error_table)
        hout.close()

if __name__ == '__main__':
    main()
