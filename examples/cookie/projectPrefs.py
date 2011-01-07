#!/usr/bin/python
# Project preferencess for making cookies

import re

def set_opt_parser(optparser):
    pass

class InitializeOptions:
    # Define which files need tweaking
    code_files = ['cookie.txt']
    file_in_suffix = ''
    file_out_suffix = ''
    label_tag = '@label@'
    pattern = re.compile('@.+?@')

    # These files are intended to make visualization with Matlab easier.
    visual_files = []
    visual_file_in_suffix = ''
    visual_file_out_suffix = ''
