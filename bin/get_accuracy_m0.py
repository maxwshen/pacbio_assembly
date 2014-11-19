# Get the overall accuracy from the blasr m0 (verbose) Blasr output between 1 large fasta file and a genome

import sys
import string
import datetime
import random
import copy
import os
import commands
import numpy as np

from collections import defaultdict

def main():
  in_f = sys.argv[1]
  scores = []
  with open(in_f) as f:
    for i, line in enumerate(f):
      if len(line) > 6 and line[6] == '%':
        scores.append(float(line.split()[-1]))
  for s in scores:
    print s
  print 'mean:', np.mean(scores), 'stdev:', np.std(scores)

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start 