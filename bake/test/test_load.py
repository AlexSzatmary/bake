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
        os.chdir('../examples/poisson/bp/')

    def test_load_file(self):
        with open('constant.bp') as f:
            self.assertEqual(
                load.load(f),
                ['@label@;constant', '@forcing_function@;1.',
                 '@ideal@;0.5 * @forcing_function@ * ((x - x_l)**2 - '
                 '(x - x_l)*(x_r - x_l)) + u_l + u_r * '
                 '(x - x_l) / (x_r - x_l)', '@n@;11', '@label@;default',
                 '@forcing_function@;1.', '@x_l@;0.', '@x_r@;1.', '@u_l@;0.',
                 '@u_r@;0', '@n@;11'])

    def tearDown(self):
        os.chdir(self.cd)

if __name__ == '__main__':
    unittest.main()
