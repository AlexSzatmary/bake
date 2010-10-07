#!/usr/bin/python

import os

pipe = os.popen("grep 'DF_Infty *=.*+.*' fit.log")
DFs = [float(line.split()[2]) for line in pipe.readlines()]
pipe.close()

pipe = os.popen("grep 'rms of residuals' fit.log")
rmses = [float(line.split()[-1]) for line in pipe.readlines()]
pipe.close()

for i in xrange(len(DFs)):
    print rmses[i]/DFs[i]
