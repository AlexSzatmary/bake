#!/usr/bin/python
# mix.py
# This has several subroutines that take bake files, do string replacements on
# them, mix and match them, and so on.

import re

class Grid:
    # This represents the space over which bake iterates and operates
#    def __init__(self, tokens, list_values, n_values, N_values, tokendict):
    def __init__(self):
        self.tokens = []
        self.list_values = []
        self.n_values = []
        self.tokendict = {}

    def replace(self, string):
        """
        Like string.replace() but replaces all tags and values for a given job

        """
        for j in range(0, len(self.tokens)):
            string = string.replace(self.tokens[j], self.values[j])
        return string

def TokenValueSubValue(values, tokendict, pattern):
    """
    This takes each token's value and expands the values by substituting in
    the values of other tokens.
    """
    for j in range(len(values)):
        foundtoken = re.search(pattern, values[j])
        # This is done iteratively so that it doesn't matter what order lines
        # appear in a bake parameter file
        while foundtoken:
            values[j] = values[j].replace(
                foundtoken.group(0), values[tokendict[foundtoken.group(0)]])
            foundtoken = re.search(pattern, values[j])
    return None


def ItRunValues(grid, pattern, slice_start=0, slice_end=0):
    """
    This iterator returns one set of values for each run; one run is given
    by each iteration of ItRunValues. These values can then be subbed in for
    the corresponding tokens.
    """
    if slice_end == 0:
        slice_end = grid.N_values
        #listi and values need to be initialized
    values = [0 for i in xrange(len(grid.tokens))]
    for list_i in ItList_i(grid.n_values, slice_start, slice_end):
        # Pick the values to be used in this run
        for j in range(len(grid.tokens)):
            values[j] = grid.list_values[j][list_i[j]]
        # Do the string replace operations on the values themselves
        TokenValueSubValue(values, grid.tokendict, pattern)
        yield values


def ItList_i(n_values, slice_start, slice_end):
    """
    This is basically a thing that encloses ItList_iNoSlice and does the
    slicing math for it.
    """
    myItList_iNoSlice = ItList_iNoSlice(n_values)
    for i in xrange(0, slice_start):
        myItList_iNoSlice.next()
    for i in xrange(slice_start, slice_end):
        yield myItList_iNoSlice.next()


def ItList_iNoSlice(n_values):
    """
    Imagine an m-dimensional grid, where m is the number of tokens; the width
    of the grid in one of these directions is the number of values for that
    token.  This iterator walks through each location in that grid; list_i
    corresponds to a single position in that grid.
    """
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


def parseBPlines(lines):
    """
    This function takes lines from a bp file, plus files from options.overwrite
    It returns a freakish tuple with:
    tokes: the list of tokens
    list_values: the list of values for a given token (these are matched by
    their indices)
    n_values is the number of values for a given index
    N_values is the repeated product of n_values; it is the total number of
    runs to do.
    tokendict gives the index of a given token
    """
    # todo convert to constructor
    # Load bp file
    grid = Grid()
    #m is the number of parameters (not the number of values for the
    #parameters)
    m = 0
    for line in lines:
        elements = line.split(';')
        if elements[0] not in grid.tokens:
            grid.tokens.append(elements[0])
            grid.list_values.append(elements[1:])
            grid.n_values.append(len(elements) - 1)
            grid.tokendict[grid.tokens[-1]] = m
            m = m + 1
        else:
            # If this key is repeated
            i = grid.tokendict[elements[0]]
            grid.list_values[i] = elements[1:]
            grid.n_values[i] = len(elements) - 1

    # Count how many runs I'm going to start
    grid.N_values = 1
    for i in grid.n_values:
        grid.N_values = grid.N_values * i
    return grid
