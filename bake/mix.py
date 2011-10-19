#!/usr/bin/python
# mix.py
# This has several subroutines that take bake files, do string replacements on
# them, mix and match them, and so on.

import re
import bakedefaults
from collections import OrderedDict


class Grid:
    """
    This represents the space over which bake iterates and operates
    """
    def __init__(self, lines):
        """
        A grid is generated from a bp file
        """
        self.table = OrderedDict()

        # Load bp file, put lines in table, where the key is the key
        # and the value is the list of values on that bp line
        for line in lines:
            elements = line.split(';')
            self.table[elements[0]] = elements[1:]

# External methods
    def mix_iterator(self):
        """
        This iterator returns one set of values for each run; one run is given
        by each iteration of grid_iterator. These values can then be subbed in
        for the corresponding keys.
        """
        self.job = OrderedDict()
        for list_i in self.grid_iterator():
            # Pick the values to be used in this run
            for (k, i) in zip(self.table.keys(), list_i):
                self.job[k] = self.table[k][i]
            # Do the string replace operations on the values themselves
            self.expand_values()
            yield self.job

    def replace(self, string):
        """
        Replaces all of the tags in string with the corresponding values for
        the current job.
        """
        # self.values is assigned in mix_iterator()
        for k, v in self.job.items():
            string = string.replace(k, v)
        return string

    def get_label(self):
        """
        Getter for the label on the current iteration.
        """
        return self.job[self.label_key]

    def set_slice(self, slice_start=0, slice_end=0):
        """
        Setter
        """
        self.slice_start = slice_start
        self.slice_end = slice_end

        # Count how many runs to start
        self.N_values = 1
        for values in self.table.values():
            self.N_values = self.N_values * len(values)

        if self.slice_end == 0:
            self.slice_end = self.N_values

    def set_key_pattern(self, key_start, key_end, f):
        """
        Setter
        """
        self.key_start = key_start
        self.key_end = key_end
        self.key_pattern = re.escape(self.key_start) + r'(.*?)' + \
            re.escape(self.key_end)
        self.label_key = (self.key_start + bakedefaults.LABEL_KEY +
                          self.key_end)
        if self.label_key not in self.table.keys():
            self.infer_label(f)

# Internal methods
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
        for k, v in self.table.items():
            # The automatic label includes any keys with multiple values
            if len(v) > 1:
                # If a key has multiple values, add both its name and its key.
                # That is, if @key1@ has multiple values, label_bits will have
                # 'key1' + '@key1@' appended. This means the label includes
                # both the key's name and the particular value it has for a
                # given job.
                label_bits.append(re.search(self.key_pattern, k).group(1) + k)
        label = '-'.join(label_bits)
        # Add the label as a key-values pair to the weird data structure
        # This is as if there were in the bp file the line,
        # label
        if not label:
            raise ValueError, "The label is blank. No label was supplied "\
                "and none can be inferred."
        self.table[self.label_key] = [label]

    def grid_iterator(self):
        """
        This is basically a thing that encloses grid_iterator_noslice and does
        the slicing math for it.
        """
        n_values = [len(v) for v in self.table.values()]
        for i in xrange(self.slice_start, self.slice_end):
            yield self.variable_base(i, n_values)

    def variable_base(self, k, base):
        """
        This gives the representation of k in the base system specified in
        `base`, where the base for each digit is different.
        For binary, base = [2, 2, ..., 2]
        For decimal, base = [10, 10, ..., 10]
        If time were represented as SS:MM:HH:DD, base = [60, 60, 24, 365]
        This order is reversed for bizarre historical reasons. This means that
        the first values in a bp file get varied more frequently.
        """
        q = []
        for b in base:
            q.append(k % b)
            k /= b
        return q

    def expand_values(self):
        """
        This takes each key's value and expands the values by substituting in
        the values of other keys.
        """
        for k, v in self.job.items():
            foundkey = re.search(self.key_pattern, v)
            # This is done iteratively so that it doesn't matter what order
            # lines appear in a bake parameter file
            while foundkey:
                v = v.replace(
                    foundkey.group(0),
                    self.job[foundkey.group(0)])
                foundkey = re.search(self.key_pattern, v)
            self.job[k] = v
