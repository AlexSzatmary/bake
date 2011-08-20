# load.py
# loads and formats bp files

import re
import os
import os.path


def load(iterator):
    """
    Loads lines from an iterator and does line parsing

    1 Handles line continuation
    2 Handles include statements
    3 Handle comments at start of line
    """
    lines = []
    for l in iterator:
	# Handle line continuation with trailing backslash
	m = re.search(r'(.*)\\\s*$', l)
	while m:
	    l = m.group(1) + iterator.next().lstrip()
	    m = re.search(r'(.*)\\\s*$', l)
	# Handle include statements
 	m = re.match(r'\s*include\(\s*([^()]+)\s*\)\s*(#.*)?$', l)
 	if m:
	    lines.extend(load_file(m.group(1)))
	    l = ''
	# Handle comments at start of line
	elif re.match(r'^\s*#', l):
	    l = ''
	if l:
	    lines.append(l.replace('\n', ''))
    return lines


def load_file(f):
    """
    Turn bp file into iterator and do load() on it.
    """
    cd = os.getcwd()
    if os.path.dirname(f):
	os.chdir(os.path.dirname(f))
    with open(os.path.basename(f)) as hin:
	lines = load(hin)
    os.chdir(cd)
    return lines
