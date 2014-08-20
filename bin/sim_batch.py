# Runs de Bruijn assembly on many simulated reads

import sys
import string
import datetime
import random
import copy
import assembly
import locAL
import os

import findKTmers
import contig_nodes
import sim_genome_getreads

from collections import defaultdict

def main():
  cnames = []
  nums = []
  for name in os.listdir('/home/jeyuan/simulation/nhoods_sim4/'):
    # print name
    if os.path.isdir('/home/jeyuan/simulation/nhoods_sim4/' + name + '/'):
      tag = name.split('_')[1]
      cnames.append(tag)
      nums.append(tag[1:])


  # ktmers = findKTmers.findKTmers('/home/mshen/research/sim/sim4/extracted_sim_reads_c10885_s1250.fasta', 24, 4)
  ktmers = findKTmers.findKTmers('/home/mshen/research/sim/sim_reads4.fasta', 24, 4)

  for i in range(len(nums)):
    print cnames[i]
    contig_nodes.test_canonical_paths('/home/jeyuan/simulation/nhoods_sim4/nhood4_' + cnames[i] + '/sim4_corr_reads_nh' + nums[i] + cnames[i] + '_siz1_i30.fasta', '/home/jeyuan/simulation/nhoods_sim4/sim_genome4_region_' + cnames[i] + '_s1250.fasta', 15, 2, 'sim_genome4_' + nums[i] + '_15.2.trim2.gv', 0, 2, ktmers)

    # diameter = contig_nodes.contig_nodes('/home/jeyuan/simulation/nhoods/nhood_' + cnames[i] + '/sim_corr_reads_nh' + nums[i] + cnames[i] + '_siz1_i30.fasta', '/home/jeyuan/simulation/nhoods/sim_genome_1m_region_' + cnames[i] + '_s1250.fasta', 15, 2, 'sim_' + nums[i] + '_15.2.trim2.gv', 0, 2)

    # sim_genome_getreads.sim_genome_getreads('/home/mshen/research/sim/sim_reads4.fasta', '/home/mshen/research/sim/sim_genome4.fasta', int(nums[i]), 1250)

    

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start