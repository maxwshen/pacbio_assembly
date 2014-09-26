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
import read_fasta

def main():
  reads = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])
  _t = int(sys.argv[4])

  gap_len = 1
  findKTgapmers(reads, _k, _t, gap_len)
  return

  # If two consecutive kt-mers are not within *width*, the region in between is uncovered
  width = 1000     # Twice of 200, which is almost the one-side width of 250 in 500bp nhood

  # genome_file = '/home/mshen/research/data/e_coli_genome.fasta'
  ktmers = findKTmers(reads, _k, _t)

  # Find blind spots
  pos, genome_len = find_genomic_positions(ktmers, genome_file, _k)
  gaps, total = find_gaps(pos, genome_len, width)

  print gaps
  print total
  return

def findKTmers(reads, _k, _t):
  # Returns a set of strings
  counts = dict()   # Key = kmer, Value = t
  h, r = read_fasta.read_fasta(reads)
  for read in r:
    for i in range(len(read) - _k + 1):
      kmer = read[i : i + _k]
      if kmer in counts:
        counts[kmer] += 1
      else:
        counts[kmer] = 1

  ans = set()
  for key, val in counts.iteritems():
    if val >= _t:
      ans.add(key)
      # print key, val
  print 'Found ' + str(len(ans)) + ' (' + str(_k) + ',' + str(_t) + ')-mers in ' + str(len(r)) + ' reads'
  return ans

def findKTgapmers(reads, _k, _t, gap_len):
  # Finds [ (k,t)-mer - gap - (k,t)-mer ]
  counts = defaultdict(list)   # Key = kmer, Value = [(read header, pos)]
  h, r = read_fasta.read_fasta(reads)
  for j in range(len(r)):
    print j
    read = r[j]
    header = h[j]
    for i in range(len(read) - _k + 1):
      kmer = read[i : i + _k]
      counts[kmer].append((header, i))

  ktgapmers = []
  for key in counts.keys():
    print key
    if len(counts[key]) >= _t:
      next = defaultdict(list)  # Key = kmer, Value = [(read header, pos)]
      for (h_name, pos) in counts[key]:
        read = r[h.index(h_name)]
        next_pos = pos + _k + gap_len
        kmer = read[next_pos : next_pos + _k]
        next[kmer].append((h_name, next_pos))
      for key2 in next.keys():
        if len(next[key2]) >= _t:
          ktgapmers.append((key, gap_len, key2))
          print len(ktgapmers)

  for k in ktgapmers:
    print k
  print len(ktgapmers)
  return

  # For each kt-mer, look at all reads (gap) positions ahead and see if kt-mer exists
  # output all gap-ktmers

def find_genomic_positions(ktmers, genome_file, _k):
  # Used to find blind spots
  h_g, r_g = read_fasta.read_fasta(genome_file)
  genome = r_g[0]
  positions = []
  for i in range(len(genome)):
    g_kmer = genome[i : i + _k]
    if g_kmer in ktmers:
      positions.append(i)
  return positions, len(genome)

def find_gaps(positions, genome_len, width):
  # Find blind spots
  prev = 0
  gaps = []   # list of tuples (start, end)
  total = 0
  for i in range(len(positions)):
    curr = positions[i]
    if curr - prev >= width:
      gaps.append((prev + width / 2, curr - width / 2))
      total += curr - prev - width
    prev = positions[i]
  return gaps, total

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start