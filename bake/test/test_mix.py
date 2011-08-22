#!/usr/bin/env python

import unittest
import load
import os
import mix


class parseBPlinesTestCase(unittest.TestCase):
    def test_overwrite(self):
        # Test that a later line overwrites a previous line
        grid = mix.parseBPlines(['@foo@;bar;baz', '@foo@;rag'])
        self.assertEqual(grid.tokens[0], '@foo@')
        self.assertEqual(grid.list_values[0][0], 'rag')


class infer_label_filename_gridTestCase(unittest.TestCase):
    def setUp(self):
        self.infer_label_grid = mix.parseBPlines(
            ['@spam@;eggs;ham', '@foo@;bar', '@baz@;rag;qux'])

    def test_infer_label(self):
        self.infer_label_grid.infer_label('corge', r'@(.*)@', '@label@')
        self.assertEqual(self.infer_label_grid.list_values[
                self.infer_label_grid.tokendict['@label@']],
                'corge-spam@spam@-baz@baz@')
        print self.infer_label_grid.list_values['@label@']


class infer_label_filename_gridTestCase(unittest.TestCase):
    def setUp(self):
        self.infer_label_grid = mix.parseBPlines(
            ['@spam@;eggs;ham', '@foo@;bar', '@baz@;rag;qux'])

    def test_infer_label(self):
        self.infer_label_grid.infer_label(None, r'@(.*)@', '@label@')
        self.assertEqual(self.infer_label_grid.list_values[
                self.infer_label_grid.tokendict['@label@']],
                'spam@spam@-baz@baz@')

if __name__ == '__main__':
    unittest.main()
