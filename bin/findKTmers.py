# findKTmers.py
#
# Finds kt-mers in a fasta file 
#
#
# /usr/bin/python2.6 findKTmers.py <fasta file> <k> <t>

import sys
import string
import datetime
import random
import copy
import os
from collections import defaultdict

def main():
  reads = sys.argv[1]
  _k = int(sys.argv[2])
  _t = int(sys.argv[3])

  kmers = findKTmers(reads, _k, _t)

def findKTmers(reads, _k, _t):
  # Returns a set of strings

  isdna = False
  counts = dict()   # Key = kmer, Value = t
  readcount = 0
  with open(reads) as f:
    for i, line in enumerate(f):
      # print i
      if isdna:
        isdna = False
        dna = line.strip()
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j+_k]
          # print kmer
          if kmer in counts:
            counts[kmer] = counts[kmer] + 1
          else:
            counts[kmer] = 1
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        isdna = True
  ans = set()

  for key, val in counts.iteritems():
    if val >= _t:
      ans.add(key)
      # print key, val
  print 'Found ' + str(len(ans)) + ' (' + str(_k) + ',' + str(_t) + ')-mers in ' + str(readcount) + ' reads'

  # print ans

  genome_file = '/home/mshen/research/data/e_coli_genome.fasta'
  with open(genome_file) as f:
    genome = ''.join(f.readlines()[1:]).translate(None, '\r\n')
  for kmer in counts.keys():
    if counts[kmer] > 3:
      print kmer, counts[kmer], genome.count(kmer)

  return ans

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start