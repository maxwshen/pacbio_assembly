# Given a single error-corrected read that hopefully corresponds to some genomic region,
# try to find reads that match this genomic region based on shared k-mers

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import read_fasta as rf

from collections import defaultdict

def main():
  ec_seq_file = sys.argv[1]
  read_file = sys.argv[2]
  _k = int(sys.argv[3])

  kmer_matching(ec_seq_file, read_file, _k)
  return

def kmer_matching(ec_seq_file, read_file, _k):
  # ec_seq_file should be a fasta file with only one sequence

  hs, rs = rf.read_fasta(ec_seq_file)
  ec_seq = rs[0]

  kmers = set()
  for i in range(len(ec_seq) - _k + 1):
    kmers.add(ec_seq[i:i + _k])

  reads = dict()    # Key = header, value = num shared kmers
  hr, rr = rf.read_fasta(read_file)
  for i in range(len(rr)):
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  for key in sorted(reads, key = reads.get, reverse = True):
    print key, reads[key]

  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start