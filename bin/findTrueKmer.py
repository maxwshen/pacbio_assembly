# findTrueKmer.py
# 
# Takes in reads. Attempts to find "true" kmers by 
# creating a graph where nodes are kmers and edges connect
# kmers that are 1 insertion or 1 deletion away from each
# other.
#
# Implementated with a table approach
#
# Genome file is used to measure accuracy.

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import GenerateIndelKmers
import locAL
import read_fasta as rf

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])
  _d = int(sys.argv[4])
  cutoff = int(sys.argv[5])     # Number of k-mers to output, ordered from highest deg
  t_atleast = int(sys.argv[6])
  if sys.argv[7] == 'True':     # Filter neighbor flag
    fn = True
  else:
    fn = False

  findTrueKmer(reads_file, genome_file, _k, _d, cutoff, t_atleast, fn)
  return 

  # Batch
  for i in range(12):
    _i = str(i)
    reads_file = '/home/mshen/research/data/high_cov/ec_reads_rh_hc_n' + _i + '_gi.fasta'
    genome_file = '/home/mshen/research/data/high_cov/ec_genome_rh_hc_n' + _i + '.fasta'
    findTrueKmer(reads_file, genome_file, _k, cutoff, t_atleast, fn)
  return

  for t_atleast in range(1, 11):
    findTrueKmer(reads_file, genome_file, _k, cutoff, t_atleast, fn)


def generateRandomKmer(_k):
  nt = ['A', 'C', 'T', 'G']
  kmer = ''
  for i in range(_k):
    r = int(np.random.randint(4, size=1))
    kmer += nt[r]
  return kmer

def graphviz(genomeKmers, degrees, kmers, cutoff):
  # genomeKmers is a set of all kmers in the reads
  # degrees is a dict. Keys = kmer, Values = degree
  # kmers is a dict. Keys = kmer, Values = t

  numToOutput = cutoff
  num = copy.copy(numToOutput)
  with open('findTrueKmers.gv', 'w') as f:
    f.write('digraph G {\n')
    for key in sorted(degrees, key=degrees.get, reverse=True):
      if num == 0:
        break
      if key in genomeKmers:
        f.write(key[:-1] + ' -> ' + key[1:] + ';\n')
      else:
        f.write(key[:-1] + ' -> ' + key[1:] + ' [color="red"];\n')
      num -= 1
    f.write('}')

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
    if len(new_degrees) > cutoff:
      break

  # print 'Filtered', count, 'kmers'
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
      # if dk != i[:-1] and dk != i[1:]:
      del_ins.add(i)
  for ik in ins_kmers:
    for i in GenerateIndelKmers.genDelKmers(ik, 1)[1]:
      # if i != ik[:-1] and i != ik[1:]:
      del_ins.add(i)

  del_ins.remove(kmer)

  neighbors = [i for i in del_ins if i in degrees]
  return neighbors

def consensus(reads):
  # Constructs a multiple alignment where each column has exactly 1 kind of nt
  # Doesn't work properly right now, need to use dynamic programming
  allnt = [list(reads[0])]
  mat = [list(reads[0])]
  for r in reads[1:]:
    print mat
    curr = []
    for i in range(len(r)):
      char = r[i]
      if char != allnt[i]:
        allnt = allnt[0:i] + [char] + allnt[i:]
        for j in range(len(mat)):
          mat[j] = mat[j][0:i] + ['-'] + mat[j][i:]  
      curr.append(char)
    mat.append(curr)

  for m in mat:
    print m

def findTrueKmer(reads_file, genome_file, _k, _d, cutoff, t_cutoff, fn):
  _krange = range(_k - _d, _k + _d + 1)
  # all_kmers : Key = kmer string, value = [t, pos1, pos2, ...]
  all_kmers = [defaultdict(list) for i in range(len(_krange))]   

  h_reads, r_reads = rf.read_fasta(reads_file)
  h_gen, r_gen = rf.read_fasta(genome_file)

  for r in r_reads:
    for j in range(len(_krange)):
      k = _krange[j]
      for h in range(len(r) - k + 1):
        kmer = r[h:h + k]
        if kmer in all_kmers[j]:
          all_kmers[j][kmer][0] += 1
          all_kmers[j][kmer].append(h)
        else:
          all_kmers[j][kmer].append(1)
          all_kmers[j][kmer].append(h)

  kmers = copy.deepcopy(all_kmers[_d])

  genome = r_gen[0]
  genomeKmers = set()
  for i in range(len(genome) - _k + 1):
    genomeKmers.add(genome[i:i+_k])

  for kmer in genomeKmers:
    if kmer not in kmers:
      kmers[kmer] = [0]

  # Find degrees of all kmers
  degrees = dict()    # Key = kmer string, Value = degree int
  for kmer in kmers:
    degree = 0
    related_kmers = [[kmer]]
    for i in range(_d):
      first = related_kmers[0]
      del_kmers = []
      for j in first:
        for k in GenerateIndelKmers.genDelKmers(j, 1)[1]:
          if k != j[:-1] and k != j[1:]:
            del_kmers.append(k)

      last = related_kmers[-1]
      ins_kmers = []
      for j in last:
        for k in GenerateIndelKmers.genInsKmers(j, 1)[1]:
          if j != k[:-1] and j != k[1:]:
            ins_kmers.append(k)
      related_kmers.append(ins_kmers)
      related_kmers.insert(0, del_kmers)

    central_t_multiplier = 2
    for i in range(len(related_kmers)):
      all_kmers_curr = all_kmers[i]
      for k_word in related_kmers[i]:
        if k_word in all_kmers_curr:
          if len(k_word) == _k:
            degree += all_kmers_curr[k_word][0] * central_t_multiplier
          else:
            degree += all_kmers_curr[k_word][0]
            kmers[kmer] += all_kmers_curr[k_word][1:]

    if kmer not in degrees:
      degrees[kmer] = degree # + kmers[kmer][0]
    else:
      degrees[kmer] += degree

  if fn:
    degrees = filter_neighbors(degrees, cutoff)

  numToOutput = cutoff
  numincorrect = 0
  current_deg = 0
  num = copy.copy(numToOutput)
  best = set()
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if kmers[key] >= t_cutoff:      # Filter by t
      if num >= 0:
        current_deg = degrees[key]
      elif degrees[key] != current_deg:
        break

      if key in genomeKmers:
        pass
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'correct'
      else:
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'incorrect'
        numincorrect += 1
      best.add(key)
      num -= 1


  # print '\nCorrect Only:'
  # for key in sorted(degrees, key=degrees.get, reverse=True):
  #   if key in genomeKmers:
  #     print key, 'Deg =', degrees[key], 't =', kmers[key], 'correct'  
  #     genomeKmers.remove(key)
  # print 'Missing:', genomeKmers

  graphviz(genomeKmers, degrees, kmers, cutoff)

  # print '\nIncorrect'
  # for key in sorted(degrees, key=degrees.get, reverse=True):
  #   if key not in genomeKmers:
  #     print key, 'Deg =', degrees[key], 't =', kmers[key]

  print 't-cutoff:', t_cutoff, 'Incorrect:', numincorrect, 'Total:', cutoff - num
  return


if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start