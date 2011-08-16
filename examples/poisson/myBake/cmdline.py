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
    # For example, the code_files string is converted to a list object
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

    # Prepends the overwrite (-o) lines from the command line to the bp file
    if options.overwrite:
        lines = options.overwrite
    else:
        lines = []

    # Scan bp file
    hin = open(options.file, 'r')
    lines += hin.readlines()
    hin.close()

    # Set up bake mixIterator
    (label, tokens, mixIterator) = bake.make_iterator(
        config['label']['label_tag'], config['label']['pattern'],
        lines, options.slice_start, options.slice_end)

    # If you were to have many custom commands, set this up like a series of
    # elif's
    if options.compare_ideal:
        # To loop over the grid of jobs, do a for loop on mixIterator
        for values in mixIterator:
            cd = values[label]
            if options.list:
                print(cd)
            wd = os.path.join('.', 'batch', cd)
            os.chdir(wd)
            h = os.popen('python compare_ideal.py')
            error = h.read()
            h.close()
            os.chdir('../..')

            table_line = '@n@ ' + error
            # This for loop is set up so that any tag could be put in
            # the default table_line and replaced for its value
            for j in range(0, len(tokens)):
                table_line = table_line.replace(tokens[j], values[j])
            error_table.append(table_line)
    # You can add other tasks like this
    # elif options.other_task:
    #     code for other task
    # bake.default_loop preserves bake's regular behavior plus this is how
    # the poisson task is run
    else:
        bake.default_loop(label, tokens, mixIterator, config, options)

    # Post-processing
    if options.compare_ideal:
        hout = open('error.txt', 'w')
        hout.writelines(error_table)
        hout.close()

if __name__ == '__main__':
    main()
