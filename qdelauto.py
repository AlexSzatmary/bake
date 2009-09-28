#!/usr/bin/python
import os, sys

for i in range(int(sys.argv[1]), int(sys.argv[2])):
    os.system('qdel ' + str(i))
