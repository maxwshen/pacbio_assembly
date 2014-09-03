# Clusters reads in a neighborhood

import sys
import string
import datetime
import random
import copy
import assembly
import os

import locAL
import read_fasta

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  
  cluster_reads(reads_file)
  return

def cluster_reads(reads_file):
  h, r = read_fasta.read_fasta(reads_file)
  mat = [[0] * len(r) for i in range(len(r))]

  for i in range(len(r)):
    mat[i][i] = 1
    for j in range(i + 1, len(r)):
      r1 = r[i]
      r2 = r[j]
      (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(r1, r2, 1, -1, -1, -0.5)
      accuracy = 0
      if alignLen > 0:
        accuracy = float(matches) / float(alignLen)
      mat[i][j] = accuracy
      mat[j][i] = accuracy
  for a in mat:
    print a

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start