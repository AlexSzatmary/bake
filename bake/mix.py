#!/usr/bin/python

# mix.py

# This has several subroutines that take bake files, do string replacements on
# them, mix and match them, and so on.
import re

# This takes each token's value and expands the values by substituting in
# the values of other tokens.
def TokenValueSubValue(values, tokendict, pattern):
    for j in range(len(values)):
        foundtoken = re.search(pattern, values[j])
        # This is done iteratively so that it doesn't matter what order lines
        # appear in a bake parameter file
        while foundtoken:
            values[j] = values[j].replace(foundtoken.group(0),
                                          values[tokendict[foundtoken.group(0)]])
            foundtoken = re.search(pattern, values[j])
    return None

# This iterator returns one set of values for each run; one run is given
# by each iteration of ItRunValues. These values can then be subbed in for
# the corresponding tokens.
def ItRunValues(list_values, tokens, n_values, N_values, pattern, tokendict,
                slice_start=0, slice_end=0):
    if slice_end == 0:
        slice_end = N_values
        #listi and values need to be initialized
    values = [0 for i in xrange(len(tokens))]
    for list_i in ItList_i(n_values, slice_start, slice_end):
        # Pick the values to be used in this run
        for j in range(len(tokens)):
            values[j] = list_values[j][list_i[j]]
        # Do the string replace operations on the values themselves
        TokenValueSubValue(values, tokendict, pattern)
        yield values

# This is basically a thing that encloses ItList_iNoSlice and does the slicing
# math for it.
def ItList_i(n_values, slice_start, slice_end):
    myItList_iNoSlice = ItList_iNoSlice(n_values)
    for i in xrange(0, slice_start):
        myItList_iNoSlice.next()
    for i in xrange(slice_start, slice_end):
        yield myItList_iNoSlice.next()

# Imagine an m-dimensional grid, where m is the number of tokens; the width of
# the grid in one of these directions is the number of values for that token.
# This iterator walks through each location in that grid; list_i corresponds
# to a single position in that grid.
def ItList_iNoSlice(n_values):
    list_i = [0 for i in xrange(len(n_values))]
    yield list_i
    while True:
        j = 0
        while True:
            list_i[j] = list_i[j] + 1
            if list_i[j] == n_values[j]:
                list_i[j] = 0
                j = j + 1
            else:
                break
        yield list_i

# This function takes lines from a bp file, plus files from options.overwrite
# It returns a freakish tuple with:
# tokes: the list of tokens
# list_values: the list of values for a given token (these are matched by their
#   indices)
# n_values is the number of values for a given index
# N_values is the repeated product of n_values; it is the total number of
#   runs to do.
# tokendict gives the index of a given token
def parseBPlines(lines):
    # Load bp file
    tokens = []
    list_values = []
    n_values = []
    tokendict = {}
    #m is the number of parameters (not the number of values for the
    #parameters)
    m = 0
    for line in lines:
        if (line[0] != '#'):
            line = line.replace('\n','').replace('\\n','&\n')
            elements = line.split(';')
            if elements[0] not in tokens:
                tokens.append(elements[0])
                list_values.append(elements[1:])
                n_values.append(len(elements)-1)
                tokendict[tokens[-1]] = m
                m = m + 1

# Count how many runs I'm going to start
    N_values = 1
    for i in n_values:
        N_values = N_values*i
    return (tokens, list_values, n_values, N_values, tokendict)
