#!/usr/bin/env python

import unittest
import load
import os

class loadLinesTestCase(unittest.TestCase):
    def test_load(self):
	# Test line continuation
	self.assertEqual(load.load(x for x in ['l1\\', '\tl2\\  ', '\tl3'])[0],
			 'l1l2l3')

class loadFileTestCase(unittest.TestCase):
    def setUp(self):
	self.cd = cd = os.getcwd()
	os.chdir('../examples/cookie/bp/')

    def test_load_file(self):
	with open('include3.txt') as f:
	    self.assertEqual(
		load.load(f),
		['@label@;include1', '@label@;cookie-decoration-@decoration@', 
		 '@type@;sugar', '@nuts@;none', '@decoration@;none', 
		 '@cutter@;circle'])

    def tearDown(self):
	os.chdir(self.cd)

if __name__ == '__main__':
    unittest.main()
