#!/usr/bin/env python

import os

bakecmd = '../../bake/cmdline.py'

p = os.popen(bakecmd + ' -l -f bp/constant.bp')
hout = open('batch/test1.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen(bakecmd + " -m -f bp/constant.bp")
hout = open('batch/test2.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen(bakecmd + " -e 'python poisson.py' -f bp/constant.bp")
hout = open('batch/test3.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen(bakecmd + " -e 'python compare_ideal.py' -f bp/constant.bp")
hout = open('batch/test4.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen(bakecmd + " -m -e 'python poisson.py; python compare_ideal.py' "
	     "-f bp/sine.bp")
hout = open('batch/test5.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen(bakecmd + " -l -m -e 'python poisson.py; "
	     "python compare_ideal.py' "
	     "-f bp/constant.bp -o '@label@;constant-n@n@' "
	     "-o '@n@;3;5;11;21;41'")
hout = open('batch/test6.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen("./myBake/cmdline.py -m -P -f bp/constant.bp "
	     "-o '@label@;mbconstant'")
hout = open('batch/test7.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

p = os.popen("./myBake/cmdline.py -l -c -f bp/constant.bp "
	     "-o '@label@;constant-n@n@' -o '@n@;3;5;11;21;41'")
hout = open('batch/test8.txt', 'w')
hout.write(p.read())
hout.close()
p.close()

# diff tests with reference case
p = os.popen("for f in `ls test`; do echo $f; diff -r batch/$f test/$f|grep -v '^Binary files'; done")
print p.read()
p.close()
