#!/usr/bin/env python

import os

bakecmd = '../../bake/cmdline.py'
dashos =  "-o '@possessive@;your' -o '@desire@;love' -o '@location@;beachhouse' -o '@adjective1@;amazing'"
my_slice = '-s 121-144'

p = os.popen(bakecmd + " -l -f trex.bp " + dashos + "|wc -l")
hout = open('batch/test-count-jobs.txt', 'w')
hout.write(p.read())
hout.close()
p.close()
p = os.popen(bakecmd + " -l -f trex.bp " + dashos)
hout = open('batch/test-list.txt', 'w')
hout.write(p.read())
hout.close()
p.close()
p = os.popen(bakecmd + " -l -f trex.bp " + my_slice + " " + dashos)
hout = open('batch/test-slice.txt', 'w')
hout.write(p.read())
hout.close()
p.close()
p = os.popen(bakecmd + " -l -f trex.bp " + my_slice + " -m -b trex.txt " + dashos)
hout = open('batch/test-mix.txt', 'w')
hout.write(p.read())
hout.close()
p.close()
p = os.popen(bakecmd + " -l -f trex.bp " + my_slice + " " + dashos + " -e 'cat trex.txt'")
hout = open('batch/test-execute.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

# diff tests with reference case
p = os.popen('for f in `ls test`; do echo $f; diff -r batch/$f test/$f; done')
print p.read()
p.close()
