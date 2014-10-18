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
  gap_len = int(sys.argv[5])
  # genome_file = '/home/mshen/research/data/e_coli_genome.fasta'

  ktgapmers = findKTgapmers(reads, _k, _t, gap_len)
  ktmers = findKTmers(reads, _k * 2, _t)

  # print ktgapmers, '\n', ktmers
  num_overlap = 0
  print 'Overlap:'
  for (k1, gap_len, k2) in ktgapmers:
    for ktmer in ktmers:
      if k1 == ktmer[:_k] and k2[:-1] == ktmer[-_k + 1:]:
        num_overlap += 1
        # print (k1, gap_len, k2), ktmer
  print num_overlap
  return

  # If two consecutive kt-mers are not within *width*, the region in between is uncovered
  width = 1000     # Twice of 200, which is almost the one-side width of 250 in 500bp nhood
  # Find blind spots
  pos, genome_len = find_genomic_positions(ktmers, genome_file, _k)
  gaps, total = find_gaps(pos, genome_len, width)

  print gaps
  print total
  return

def findKTmers(reads, _k, _t):
  # Returns a set of strings
  isdna = False
  counts = dict()
  r = 0
  with open(reads) as f:
    for i, line in enumerate(f):
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
        r += 1
        isdna = True

  ans = set()
  for key, val in counts.iteritems():
    if val >= _t:
      ans.add(key)
      # print key, val
  print 'Found ' + str(len(ans)) + ' (' + str(_k) + ',' + str(_t) + ')-mers in ' + str(r) + ' reads'
  return ans

def findKTgapmers(reads, _k, _t, gap_len):
  # Finds [ (k,t)-mer - gap - (k,t)-mer ]
  isdna = False
  counts = defaultdict(list)
  r = 0
  with open(reads) as f:
    for i, line in enumerate(f):
      if isdna:
        isdna = False
        dna = line.strip()
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j + _k]
          kmer2 = dna[j + _k + gap_len : j + 2 * _k + gap_len]
          if len(kmer2) == _k:
            counts[kmer].append(kmer2)
      if line[0] == '>' or line[0] == '@':
        r += 1
        isdna = True


  ktgapmers = set()
  for kmer in counts.keys():
    for kmer2 in counts[kmer]:
      if counts[kmer].count(kmer2) >= _t:
        ktgapmers.add((kmer, gap_len, kmer2))

  # for key in counts.keys():
  #   if len(counts[key]) >= _t:
  #     next = defaultdict(list)  # Key = kmer, Value = [(read header, pos)]
  #     for (h_name, pos) in counts[key]:
  #       read = find_read(reads, h_name)
  #       next_pos = pos + _k + gap_len
  #       kmer = read[next_pos : next_pos + _k]
  #       if len(kmer) == _k:
  #         next[kmer].append((h_name, next_pos))
  #     # print key, next
  #     for key2 in next.keys():
  #       if len(next[key2]) >= _t:
  #         ktgapmers.append((key, gap_len, key2))

  print ktgapmers
  print 'Found ' + str(len(ktgapmers)) + ' (' + str(_k) + ',' + str(_t) + ',' + str(gap_len) + ')-mers in ' + str(r) + ' reads'
  return ktgapmers

  # For each kt-mer, look at all reads (gap) positions ahead and see if kt-mer exists
  # output all gap-ktmers

def find_read(reads_file, header):
  found = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if found:
        return line.strip()
      if line == header:
        found = True

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