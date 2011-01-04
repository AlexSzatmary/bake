#!/usr/bin/python

import os
try:
    p = os.popen('../../bake.py -m bp/many_sprinkles.txt')
    hout = open('batch/test1.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake.py -l bp/many_sprinkles.txt -s 2-4')
    hout = open('batch/test2.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake.py -m bp/many_sprinkles.txt -o "@color@;purple"')
    hout = open('batch/test3.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('for f in `ls batch`; do echo $f; diff -r batch/$f test/$f; done')
    print p.read()
    p.close()
except:
    pass
