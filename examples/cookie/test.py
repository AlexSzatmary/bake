#!/usr/bin/python

import os
try:
    p = os.popen('../../bake.py -m -l -f bp/many_sprinkles.txt')
    hout = open('batch/test1.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake.py -l -f bp/many_sprinkles.txt -s 2-4')
    hout = open('batch/test2.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake.py -m -l -f bp/many_sprinkles.txt -o "@color@;purple"')
    hout = open('batch/test3.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close() 
    p = os.popen('../../bake.py -f bp/many_sprinkles.txt -e "cat cookie.txt"')
    hout = open('batch/test4.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    

    # diff tests with reference case
    p = os.popen('for f in `ls batch`; do echo $f; diff -r batch/$f test/$f; done')
    print p.read()
    p.close()
except:
    pass
