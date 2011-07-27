#!/usr/bin/python

import re

def TokenValueSubValue(values, tokendict, pattern):
  for j in range(len(values)):
    foundtoken = re.search(pattern, values[j])
    while foundtoken:
      values[j] = values[j].replace(foundtoken.group(0), 
                                    values[tokendict[foundtoken.group(0)]])
      foundtoken = re.search(pattern, values[j])
  return None

def ItRunValues(list_values, tokens, n_values, N_values, pattern, tokendict,
                slice_start=0, slice_end=0):
  if slice_end == 0:
    slice_end = N_values
    #listi and values need to be initialized
  values = [0 for i in xrange(len(tokens))]
  for list_i in ItList_i(n_values, slice_start, slice_end):
    # Pick the values to be used in this run
    for j in range(len(tokens)):
      values[j] = list_values[j][list_i[j]]
    # Do the string replace operations on the values themselves
    TokenValueSubValue(values, tokendict, pattern)
    yield values

def ItList_i(n_values, slice_start, slice_end):
  myItList_iNoSlice = ItList_iNoSlice(n_values)
  for i in xrange(0, slice_start):
    myItList_iNoSlice.next()
  for i in xrange(slice_start, slice_end):
    yield myItList_iNoSlice.next()

def ItList_iNoSlice(n_values):
    list_i = [0 for i in xrange(len(n_values))]
    yield list_i
    while True:
      j = 0
      while True:
	list_i[j] = list_i[j] + 1
	if list_i[j] == n_values[j]:
	  list_i[j] = 0
	  j = j + 1
	else:
	  break
      yield list_i

def parseBPlines(lines):
    # Load bp file
    tokens = []
    list_values = []
    n_values = []
    tokendict = {}
    #m is the number of parameters (not the number of values for the 
    #parameters)
    m = 0
    for line in lines:
      if (line[0] != '#'):
	line = line.replace('\n','').replace('\\n','&\n')
	elements = line.split(';')
        if elements[0] not in tokens:
          tokens.append(elements[0])
          list_values.append(elements[1:])
          n_values.append(len(elements)-1)
          tokendict[tokens[-1]] = m
          m = m + 1

# Count how many runs I'm going to start
    N_values = 1
    for i in n_values:
      N_values = N_values*i
    return (tokens, list_values, n_values, N_values, tokendict)
