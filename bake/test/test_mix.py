#!/usr/bin/env python

import unittest
import load
import os
import mix


class parseBPlinesTestCase(unittest.TestCase):
    def test_overwrite(self):
        # Test that a later line overwrites a previous line
        grid = mix.Grid(['@foo@;bar;baz', '@foo@;rag'])
        self.assertEqual(grid.table['@foo@'], ['rag'])


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
        self.infer_label_grid.set_key_pattern('@', '@', None)
        self.assertEqual(self.infer_label_grid.table['@label@'],
                         ['spam@spam@-baz@baz@'])

if __name__ == '__main__':
    unittest.main()
