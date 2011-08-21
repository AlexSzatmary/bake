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

if __name__ == '__main__':
    unittest.main()
