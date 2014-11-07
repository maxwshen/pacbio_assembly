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
from subprocess import call

def main():
  # reads_file = sys.argv[1]
  directory = sys.argv[1]

  batch(directory)

  # cluster_reads(reads_file)
  return

def batch(directory):
  files = os.listdir(directory)
  for fil in files:
    print directory + fil, datetime.datetime.now()
    cluster_reads(directory + fil)

def cluster_reads(reads_file):
  h, r = read_fasta.read_fasta(reads_file)
  mat = [[0] * len(r) for i in range(len(r))]

  len_cutoff = 400
  r = filter_by_len(r, len_cutoff)

  clusters = []

  while len(r) > 0:
    r1 = r[0]
    curr_cluster = [h[0], r1]
    for j in range(1, len(r)):
      r2 = r[j]
      (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external(r1, r2, 1, -2, -2, -1)
      # print j, alignLen
      if alignLen > len_cutoff:
        curr_cluster.append(h[j])
        curr_cluster.append(r[j])
    clusters.append(curr_cluster)
    for item in curr_cluster:
      if item[0] == '>':
        del r[h.index(item)]
        del h[h.index(item)]
    # print curr_cluster, '\n', len(curr_cluster)
    # print len(r), len(h)    

  for i in range(len(clusters)):
    out_file = 'clusters_22.4_100k/nhood_' + reads_file.split('_')[6] + '_cluster_' + str(i) + '.fasta'
    with open(out_file, 'w') as f:
      for line in clusters[i]:
        f.write(line + '\n')

  return


def filter_by_len(reads, cutoff):
  new_r = []
  for r in reads:
    if len(r) > cutoff:
      new_r.append(r)
  return new_r

def longest_common_substring(s1, s2):
   m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
   longest, x_longest = 0, 0
   for x in xrange(1, 1 + len(s1)):
       for y in xrange(1, 1 + len(s2)):
           if s1[x - 1] == s2[y - 1]:
               m[x][y] = m[x - 1][y - 1] + 1
               if m[x][y] > longest:
                   longest = m[x][y]
                   x_longest = x
           else:
               m[x][y] = 0
   return s1[x_longest - longest: x_longest]

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start