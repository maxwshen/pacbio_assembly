# compare_nhood.py
# Given a neighborhood, finds the corresponding true reads in the same region
# and compares the two sets of reads

import sys
import string
import datetime
import random
import copy
import os
from collections import defaultdict

import read_fasta
import PacBioRegionStats as pbrs

def main():
  nhood_file = sys.argv[1]
  _k = int(sys.argv[2])
  _t = int(sys.argv[3])

  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_whole_removed_homopolymers.fasta'
  genome_file = '/home/mshen/research/data/e_coli_genome.fasta'

  compare_nhood(nhood_file, reads_file, genome_file, _k, _t)

def compare_nhood(nhood_file, reads_file, genome_file, _k, _t):
  # Expected nhood_file format: nhood_nh12345_ACTGACTG.fasta
  #
  # Sample read header:
  # >m120114_011938_42177_c100247042550000001523002504251220_s1_p0/3500/0_6472/0_6472
  
  pos = int(nhood_file.split('/')[-1].split('_')[1][2:])
  width = 500
  min_overlap = 150

  extract_reads_f, extract_genome_f = pbrs.extract(pos, width, min_overlap, reads_file, genome_file)

  nhood_h, nhood_r = read_fasta.read_fasta(nhood_file)
  reads_h, reads_r = read_fasta.read_fasta(extract_reads_f)

  reads_h_trim = []
  for h in reads_h:
    reads_h_trim.append(h[1:].split('/')[:2])

  in_reads = []
  not_in_reads = []
  for h in nhood_h:
    h_trim = h[1:].split('/')[:2]
    if h_trim in reads_h_trim:
      in_reads.append(h_trim)
    else:
      not_in_reads.append(h_trim)

  # print in_reads, not_in_reads
  graph_ktmers(nhood_h, nhood_r, _k, _t)
  print 'In reads:', len(in_reads), '\tNot:', len(not_in_reads), '\tMissing:', len(reads_h) - len(in_reads)
  print '|Nhood|:', len(nhood_h), '\t|Reads|:', len(reads_h)
  return

def graph_ktmers(nhood_h, nhood_seqs, _k, _t):
  ktmers = dict()   # Key = kmer, Value = multiplicity
  for seq in nhood_seqs:
    for i in range(len(seq) - _k + 1):
      kmer = seq[i : i + _k]
      if kmer in ktmers:
        ktmers[kmer] += 1
      else:
        ktmers[kmer] = 1
  nodes = [k for k in ktmers.keys() if ktmers[k] >= _t]

  all_vis = ''
  sum_ktmer_occur = 0
  for j in range(len(nhood_seqs)):
    seq = nhood_seqs[j]
    header = nhood_h[j]
    vis = ''
    prev = 0
    for i in range(len(seq)):
      kmer = seq[i : i + _k]
      if kmer in nodes:
        sum_ktmer_occur += 1
        vis += str(i - prev) + ' - ' + kmer + ' - '
        prev = i
    vis += str(len(seq) - prev)
    all_vis += header + '\n' + vis + '\n'
  print all_vis
  print 'Avg. # ktmers/read:', float(sum_ktmer_occur) / float(len(nhood_seqs))
  return


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start