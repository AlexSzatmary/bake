#!/usr/bin/python

import time, os, random

random.seed()
runs = 'foo'
tnorm = 60*60

while runs != '':
    if runs != 'foo':
	time.sleep(random.uniform(tnorm*0.75, tnorm*1.25))
    runsssh = os.popen(
        'ssh al1@hpc.rs.umbc.edu qstat|grep al1')
    runs = runsssh.read()
#    print runs
    runsssh.close()

os.system("growlnotify -n Runmonitor -s -m 'hpc runs done'")
