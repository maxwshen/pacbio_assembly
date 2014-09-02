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
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])

  highdegnodes_stats(reads_file, genome_file, _k)

def highdegnodes_stats(reads_file, genome_file, _k):
  h_reads, r_reads = read_fasta(reads_file)
  h_genome, r_genome = read_fasta(genome_file)
  genome = r_genome[0]

  Kmers = []

  for x in range(len(genome) - _k + 1):
    kmer = genome[x : x + _k]
    Kmers.append(Kmer(kmer))
    curr = Kmers[-1]
    curr.pos_true = x

    # Remove ends
    del_kmers = GenerateIndelKmers.genDelKmers(kmer, 1)[1]
    remove_list = [dk for dk in del_kmers if dk == kmer[:-1] or dk == kmer[1:]]
    for dk in remove_list:
      del_kmers.remove(dk)
    ins_kmers = GenerateIndelKmers.genInsKmers(kmer, 1)[1]
    remove_list = [ik for ik in ins_kmers if kmer == ik[:-1] or kmer == ik[1:]]
    for ik in remove_list:
      ins_kmers.remove(ik)

    # Deletion and insertion
    del_ins = set()
    for dk in del_kmers:
      for i in GenerateIndelKmers.genInsKmers(dk, 1)[1]:
        del_ins.add(i)
    for ik in ins_kmers:
      for i in GenerateIndelKmers.genDelKmers(ik, 1)[1]:
        del_ins.add(i)
    del_ins.remove(kmer)

    ins_ins = set()
    del_del = set()
    for ik in ins_kmers:
      for i in GenerateIndelKmers.genInsKmers(ik, 1)[1]:
        ins_ins.add(i)
    for dk in del_kmers:
      for i in GenerateIndelKmers.genDelKmers(dk, 1)[1]:
        del_del.add(i)

    print '\n', kmer, x

    for j in range(len(r_reads)):
      read = r_reads[j]
      starting_pos = int(h_reads[j].split('$')[-3])
      # Format: >m120114_011938_42177_c100247042550000001523002504251220_s1_p0/756/0_7141/692_1792/$0$971$/

      for i in range(len(read)):
        l_0 = read[i : i + _k]
        l_p1 = read[i : i + _k + 1]
        l_p2 = read[i : i + _k + 2]
        l_m1 = read[i : i + _k - 1]
        l_m2 = read[i : i + _k - 2]

        pos = starting_pos + i

        # Exact
        if l_0 == kmer:
          curr._t += 1
        if l_0 in del_ins:
          curr.pos_0.append(pos)
          print l_0, pos
        if l_p1 in ins_kmers:
          curr.pos_p1.append(pos)
        if l_m1 in del_kmers:
          curr.pos_m1.append(pos)
        if l_p2 in ins_ins:
          curr.pos_p2.append(pos)
        if l_m2 in del_del:
          curr.pos_m2.append(pos)
  
  print 'kmer\tk+1\tk+2\tk\tk-1\tk-2\tavg pos\ttrue pos\tdiff\tmax.pos\tmin.pos'
  for k in Kmers:
    print k

class Kmer():
  def __init__(self, kmer):
    self.kmer = kmer
    self.pos_p1 = []
    self.pos_p2 = []
    self.pos_0 = []
    self.pos_m1 = []
    self.pos_m2 = []
    self.pos_true = -1
    self._t = 1

  def __str__(self):
    all_pos = [] + self.pos_p1 + self.pos_p2 + self.pos_0 + self.pos_m1 + self.pos_m2
    return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}'.format(self.kmer, len(self.pos_p1), len(self.pos_p2), len(self.pos_0), len(self.pos_m1), len(self.pos_m2), np.mean(all_pos), self.pos_true, np.absolute(np.mean(all_pos) - self.pos_true), max(all_pos), min(all_pos))



def read_fasta(fasta_file):
  # Builds header and reads lists from a fasta file
  headers = []
  reads = []
  get_read = False
  with open(fasta_file) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] == '>' or line[0] == '@':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        headers.append(line.strip())
        get_read = True
      elif get_read:
        curr_read += line.strip()
  if curr_read != '':
    reads.append(curr_read)
  if len(headers) == 0 or len(reads) == 0:
    print 'ERROR: Bad fasta file', fasta_file
    sys.exit(0)
  return headers, reads

if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start