# findHighDegreeKmers.py
# 
# reads_file: Fasta format reads
# k: The length of the k-mer to find
# cutoff: How many k-mers to output, ordered from highest degree to lowest
#
# Returns a set of k-mers

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import GenerateIndelKmers
import locAL

from collections import defaultdict

def main():
  if len(sys.argv) != 4:
    print 'Usage: python findHighDegreeKmers <reads_file> <k> <cutoff>'
    sys.exit(0)

  findHighDegreeKmers(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))

def generateRandomKmer(_k):
  nt = ['A', 'C', 'T', 'G']
  kmer = ''
  for i in range(_k):
    r = int(np.random.randint(4, size=1))
    kmer += nt[r]
  return kmer

def findHighDegreeKmers(reads_file, _k, cutoff):
  _kplus = _k + 1
  _kminus = _k - 1
  isdna = False
  readcount = 0
  kmers = dict()
  kplusmers = dict()
  kminusmers = dict()

  with open(reads_file) as f:
    for i, line in enumerate(f):
      if isdna:
        isdna = False
        dna = line.strip()
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j+_k]
          if kmer in kmers:
            kmers[kmer] += 1
          else:
            kmers[kmer] = 1
        for j in range(len(dna) - _kminus + 1):
          kmer = dna[j:j+_kminus]
          if kmer in kminusmers:
            kminusmers[kmer] += 1
          else:
            kminusmers[kmer] = 1
        for j in range(len(dna) - _kplus + 1):
          kmer = dna[j:j+_kplus]
          if kmer in kplusmers:
            kplusmers[kmer] += 1
          else:
            kplusmers[kmer] = 1
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        isdna = True

  degrees = dict()
  for kmer in kmers:
    del_kmers = GenerateIndelKmers.genDelKmers(kmer, 1)[1]
    degree = 0
    for del_kmer in del_kmers:
      if del_kmer in kminusmers:
        # degree += 1
        degree += kminusmers[del_kmer]

    ins_kmers = GenerateIndelKmers.genInsKmers(kmer, 1)[1]
    for ins_kmer in ins_kmers:
      if ins_kmer in kplusmers:
        # degree += 1 
        degree += kplusmers[ins_kmer]

    degrees[kmer] = degree

  numToOutput = cutoff
  num = copy.copy(numToOutput)
  best = set()
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if num == 0:
      break
    print key, 'Deg =', degrees[key], 't =', kmers[key]
    best.add(key)
    num -= 1

  return best


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start