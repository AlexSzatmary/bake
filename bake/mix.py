#!/usr/bin/python
# mix.py
# This has several subroutines that take bake files, do string replacements on
# them, mix and match them, and so on.

import re
import bakedefaults


class Grid:
    """
    This represents the space over which bake iterates and operates
    """
    def __init__(self, lines):
        """
        A grid is generated from a bp file

        It has:
        keys: the list of keys
        list_values: the list of values for a given key (these are matched by
        their indices)
        N_values is the repeated product of n_values; it is the total number of
        runs to do.
        keydict gives the index of a given key
        """
        self.keys = []
        self.list_values = []
        self.keydict = {}

        # Load bp file
        # m is the number of keys (not the number of values for the
        # keys)
        m = 0
        for line in lines:
            elements = line.split(';')
            if elements[0] not in self.keys:
                self.keys.append(elements[0])
                self.list_values.append(elements[1:])
                self.keydict[self.keys[-1]] = m
                m = m + 1
            else:
                # If this key is repeated
                i = self.keydict[elements[0]]
                self.list_values[i] = elements[1:]

        # Count how many runs I'm going to start
        self.N_values = 1
        for values in self.list_values:
            self.N_values = self.N_values * len(values)

    def infer_label(self, string):
        """
        Makes up a sensible label if none is specified.

        string should be something like the bp filename, a good start to the
        label. If it's absent, that doesn't matter.

        The label format is:
        string-key1@key1@-key2@key2...
        or
        key1@key1@-key2@key2...
        if string is absent.
        """
        label_bits = []
        if string:
            label_bits.append(string)
        for i in xrange(len(self.list_values)):
            # The automatic label includes any keys with multiple values
            if len(self.list_values[i]) > 1:
                # If a key has multiple values, add both its name and its key.
                # That is, if @key1@ has multiple values, label_bits will have
                # 'key1' + '@key1@' appended. This means the label includes
                # both the key's name and the particular value it has for a
                # given job.
                label_bits.append(
                    re.search(self.key_pattern, self.keys[i]).group(1) +
                    self.keys[i])
        label = '-'.join(label_bits)
        # Add the label as a key-values pair to the weird data structure
        # This is as if there were in the bp file the line,
        # label
        if not label:
            raise ValueError, "The label is blank. No label was supplied "\
                "and none can be inferred."
        self.keys.append(self.label_key)
        self.list_values.append([label])
        self.keydict[self.label_key] = len(self.list_values) - 1

    def replace(self, string):
        """
        Replaces all of the tags in string with the corresponding values for
        the current job.
        """
        # self.values is assigned in mix_iterator()
        for j in xrange(len(self.list_values)):
            string = string.replace(self.keys[j], self.values[j])
        return string

    def expand_values(self):
        """
        This takes each key's value and expands the values by substituting in
        the values of other keys.
        """
        for j in range(len(self.values)):
            foundkey = re.search(self.key_pattern, self.values[j])
            # This is done iteratively so that it doesn't matter what order
            # lines appear in a bake parameter file
            while foundkey:
                self.values[j] = self.values[j].replace(
                    foundkey.group(0),
                    self.values[self.keydict[foundkey.group(0)]])
                foundkey = re.search(self.key_pattern, self.values[j])

    def set_slice(self, slice_start=0, slice_end=0):
        """
        Setter
        """
        self.slice_start = slice_start
        self.slice_end = slice_end
        if self.slice_end == 0:
            self.slice_end = self.N_values

    def set_key_pattern(self, args):
        self.key_start = args.key_start
        self.key_end = args.key_end
        self.key_pattern = re.escape(self.key_start) + r'(.*?)' + \
            re.escape(self.key_end)
        self.label_key = (self.key_start + bakedefaults.LABEL_KEY +
                          self.key_end)
        if self.label_key not in self.keys:
            self.infer_label(args.file)

    def mix_iterator(self):
        """
        This iterator returns one set of values for each run; one run is given
        by each iteration of ItRunValues. These values can then be subbed in
        for the corresponding keys.
        """
            #listi and values need to be initialized
        self.values = [0 for i in xrange(len(self.keys))]
        for list_i in self.grid_iterator():
            # Pick the values to be used in this run
            for j in range(len(self.keys)):
                self.values[j] = self.list_values[j][list_i[j]]
            # Do the string replace operations on the values themselves
            self.expand_values()
            yield self.values

    def get_label(self):
        """
        Getter for the label on the current iteration.
        """
        return self.values[self.keydict[self.label_key]]

    def grid_iterator(self):
        """
        This is basically a thing that encloses ItList_iNoSlice and does the
        slicing math for it.
        """
        my_grid_iterator_noslice = self.grid_iterator_noslice()
        for i in xrange(0, self.slice_start):
            my_grid_iterator_noslice.next()
        for i in xrange(self.slice_start, self.slice_end):
            yield my_grid_iterator_noslice.next()

    def grid_iterator_noslice(self):
        """
        Imagine an m-dimensional grid, where m is the number of keys; the width
        of the grid in one of these directions is the number of values for that
        key. This iterator walks through each location in that grid; list_i
        corresponds to a single position in that grid.
        """
        list_i = [0 for i in xrange(len(self.list_values))]
        yield list_i
        while True:
            j = 0
            while True:
                list_i[j] = list_i[j] + 1
                if list_i[j] == len(self.list_values[j]):
                    list_i[j] = 0
                    j = j + 1
                else:
                    break
            yield list_i
