#!/usr/bin/env python

import unittest
import load
import os
import mix


class parseBPlinesTestCase(unittest.TestCase):
    def test_overwrite(self):
        # Test that a later line overwrites a previous line
        grid = mix.Grid(['@foo@;bar;baz', '@foo@;rag'])
        self.assertEqual(grid.keys[0], '@foo@')
        self.assertEqual(grid.list_values[0][0], 'rag')


class infer_label_filename_gridTestCase(unittest.TestCase):
    def setUp(self):
        self.infer_label_grid = mix.Grid(
            ['@spam@;eggs;ham', '@foo@;bar', '@baz@;rag;qux'])

    def test_infer_label(self):
        self.infer_label_grid.infer_label('corge', r'@(.*)@', '@label@')
        self.assertEqual(self.infer_label_grid.get_label(),
                'corge-spam@spam@-baz@baz@')


class infer_label_filename_gridTestCase(unittest.TestCase):
    def setUp(self):
        self.infer_label_grid = mix.Grid(
            ['@spam@;eggs;ham', '@foo@;bar', '@baz@;rag;qux'])

    def test_infer_label(self):
        self.infer_label_grid.set_key_pattern({}, None)
        self.infer_label_grid.infer_label(None)
        self.assertEqual(self.infer_label_grid.list_values[
                self.infer_label_grid.keydict[
                self.infer_label_grid.label_key]][0],
                         'spam@spam@-baz@baz@')

if __name__ == '__main__':
    unittest.main()
