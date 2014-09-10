import sys
import string
import datetime
import random
import copy
import os
import numpy as np

def main():
  reads_file = sys.argv[1]
  reads_stats(reads_file)
  return

def reads_stats(reads_file):
  lengths = []
  isDNA = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if isDNA:
        lengths.append(len(line))
        print len(line)
        isDNA = False        
      if line[0] == '@':
        isDNA = True

  print 'file:', reads_file
  print 'avg len:', np.mean(lengths)
  print 'std len:', np.std(lengths)

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start