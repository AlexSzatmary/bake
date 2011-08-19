#!/usr/bin/python

import os
try:
    p = os.popen('../../bake/cmdline.py -m -l -f bp/many_sprinkles.txt')
    hout = open('batch/test1.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake/cmdline.py -l -f bp/many_sprinkles.txt -s 1-3')
    hout = open('batch/test2.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake/cmdline.py -m -l -f bp/many_sprinkles.txt -o "@color@;purple"')
    hout = open('batch/test3.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close() 
    p = os.popen('../../bake/cmdline.py -f bp/many_sprinkles.txt -e "cat cookie.txt"')
    hout = open('batch/test4.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    p = os.popen('../../bake/cmdline.py -m -l -f bp/nuts.txt -e "cat cookie.txt"')
    hout = open('batch/test5.txt', 'w')
    hout.write(p.read())
    hout.close()
    p.close()
    
    # diff tests with reference case
    p = os.popen('for f in `ls test`; do echo $f; diff -r batch/$f test/$f; done')
    print p.read()
    p.close()
except:
    pass
