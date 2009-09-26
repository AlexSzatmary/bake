#!/usr/bin/python

import time, os, random

random.seed()
runs = 'foo'
tnorm = 60*60

while runs != '':
    if runs != 'foo':
	time.sleep(random.uniform(tnorm*0.75, tnorm*1.25))
    runsssh = os.popen(
        'ssh szatmary@pople.psc.edu qstat|grep szatmary')
    runs = runsssh.read()
    print runs
    runsssh.close()

os.system("growlnotify -n Runmonitor -s -m 'Pople runs done'")
