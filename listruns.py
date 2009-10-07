#!/usr/bin/python

import re

def TokenValueSubValue(values, tokendict, pattern):
  for j in range(len(values)):
    foundtoken = re.search(pattern, values[j])
    while foundtoken:
      print foundtoken.group(0)
      print tokendict[foundtoken.group(0)]
      print values[tokendict[foundtoken.group(0)]]
      values[j] = values[j].replace(foundtoken.group(0), 
                                    values[tokendict[foundtoken.group(0)]])
      foundtoken = re.search(pattern, values[j])
    return None
