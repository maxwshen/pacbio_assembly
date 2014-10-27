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
import read_fasta

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])
  _d = int(sys.argv[4])                   # 1 or 2
  cutoff = int(sys.argv[5])               # Set as -1 to print all
  check_correctness_flag = sys.argv[6]

  findHighDegreeKmers(reads_file, genome_file, _k, _d, cutoff, check_correctness = check_correctness_flag)

def generateRandomKmer(_k):
  nt = ['A', 'C', 'T', 'G']
  kmer = ''
  for i in range(_k):
    r = int(np.random.randint(4, size=1))
    kmer += nt[r]
  return kmer

def findHighDegreeKmers(reads_file, genome_file, _k, _d, cutoff, check_correctness = False):
  # genome_file = '/home/mshen/research/extracts/extracted_genome_c3005150_s1000.fasta'
  # genome_file = '/home/mshen/research/extracts_100k/extracted_genome_c2500000_s100000.fasta'
  # genome_file = '/home/mshen/research/data/high_cov/ec_genome_rh_hc_n0.fasta'

  _kplus = _k + 1
  _kminus = _k - 1
  _kminus2 = _k - 2
  _kplus2 = _k + 2
  isdna = False
  readcount = 0
  kmers = dict()
  kplusmers = dict()
  kminusmers = dict()
  kminus2mers = dict()
  kplus2mers = dict()

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
        if _d == 2:
          for j in range(len(dna) - _kplus2 + 1):
            kmer = dna[j:j+_kplus2]
            if kmer in kplus2mers:
              kplus2mers[kmer] += 1
            else:
              kplus2mers[kmer] = 1
          for j in range(len(dna) - _kminus2 + 1):
            kmer = dna[j:j+_kminus2]
            if kmer in kminus2mers:
              kminus2mers[kmer] += 1
            else:
              kminus2mers[kmer] = 1
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        isdna = True

  degrees = dict()
  counter = 0
  for kmer in kmers:
    counter += 1
    # print counter, len(kmers)
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
    if _d == 2:
      for dk in del_kmers:
        del2_kmers = GenerateIndelKmers.genDelKmers(dk, 1)[1]
        for dk2 in del2_kmers:
          if dk2 in kminus2mers:
            degree += kminus2mers[dk2]
      for ik in ins_kmers:
        ins2_kmers = GenerateIndelKmers.genInsKmers(ik, 1)[1]
        for ik2 in ins2_kmers:
          if ik2 in kplus2mers:
            degree += kplus2mers[ik2]
    degrees[kmer] = degree + kmers[kmer]

  if check_correctness:
    h, r = read_fasta.read_fasta(genome_file)
    genome = r[0]
    genome_kmers = set()
    for i in range(0, len(genome) - _k + 1):
      genome_kmers.add(genome[i:i + _k])



  if cutoff == -1:
    numToOutput = len(degrees)
  else:
    numToOutput = cutoff
  num = copy.copy(numToOutput)

  fn = True
  if fn:
    degrees = filter_neighbors(degrees, numToOutput)

  best = set()
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if num == 0:
      break
    if check_correctness:
      if key in genome_kmers:
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'True'
      else:
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'False'
    else:
        print key, 'Deg =', degrees[key], 't =', kmers[key]
    best.add(key)
    num -= 1

  return best

def filter_neighbors(degrees, cutoff):
  # Input: Dict, keys = kmers, values = degrees
  new_degrees = dict()
  removed = set()
  count = 0
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if key not in removed:
      new_degrees[key] = degrees[key]
      neighbors = find_neighbors(degrees, key)
      for n in [i for i in neighbors if i in degrees]:
        removed.add(n)
    else:
      count += 1
    if len(new_degrees) >= cutoff:
      break

  print 'Filtered', count, 'kmers'
  return new_degrees

def find_neighbors(degrees, kmer):
  del_kmers = GenerateIndelKmers.genDelKmers(kmer, 1)[1]
  remove_list = []
  for dk in del_kmers:
    if dk == kmer[:-1] or dk == kmer[1:]:
      remove_list.append(dk)
  for dk in remove_list:
    del_kmers.remove(dk)
  remove_list = []
  ins_kmers = GenerateIndelKmers.genInsKmers(kmer, 1)[1]
  for ik in ins_kmers:
    if kmer == ik[:-1] or kmer == ik[1:]:
      remove_list.append(ik)
  for ik in remove_list:
    ins_kmers.remove(ik)

  del_ins = set()
  for dk in del_kmers:
    for i in GenerateIndelKmers.genInsKmers(dk, 1)[1]:
      del_ins.add(i)
      # if dk != i[:-1] and dk != i[1:]:
        # del_ins.add(i)

  for ik in ins_kmers:
    for i in GenerateIndelKmers.genDelKmers(ik, 1)[1]:
      del_ins.add(i)  
      # if i != ik[:-1] and i != ik[1:]:
        # del_ins.add(i)

  # del_ins.remove(kmer)

  neighbors = [i for i in del_ins if i in degrees]
  # if len(neighbors) > 2:
  #   consensus = sma.spec_multi_align(neighbors)
  #   print 'actual kmer:', kmer, consensus == kmer
  # else:
  #   for n in neighbors:
  #     print n
  #   print 'insufficient num of neighbors'
  return neighbors


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start