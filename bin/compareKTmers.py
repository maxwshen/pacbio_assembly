# Expected input is the output from findKTmers.py
# Format on each line:
# <ktmer> <T>

import sys
import string
import datetime
import random
import copy
import os
from collections import defaultdict

def main():
  if len(sys.argv) != 3:
    print 'Usage: python compareKTmers <ktmers1> <ktmers2>'
    sys.exit(0)

  compareKTmers(sys.argv[1], sys.argv[2])


def compareKTmers(ktmers1, ktmers2):
  ktmers = set()
  with open(ktmers1) as f:
    for i, line in enumerate(f):
      words = line.split()
      kmer = words[0].translate(None, ',')
      if kmer != 'Found':
        ktmers.add(kmer)

  overlap = 0
  onlyInFile2 = 0
  ktmersInFile2 = set()
  with open(ktmers2) as f:
    for i, line in enumerate(f):
      words = line.split()
      kmer = words[0].translate(None, ',')
      if kmer != 'Found':
        if kmer in ktmers:
          overlap += 1
        else:
          onlyInFile2 += 1
          ktmersInFile2.add(kmer)

  onlyInFile1 = len(ktmers) - overlap

  print 'Overlap:', overlap
  print 'Only in', ktmers1, ':', onlyInFile1
  print 'Only in', ktmers2, ':', onlyInFile2

  print 'Kmers only in', ktmers2, ':'
  for kmer in ktmersInFile2:
    print ' ', kmer

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start