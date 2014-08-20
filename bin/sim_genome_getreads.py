# Runs de Bruijn assembly on many simulated reads

import sys
import string
import datetime
import random
import copy
import assembly
import locAL
import os
import contig_nodes

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  target = int(sys.argv[3])
  width = int(sys.argv[4])

  sim_genome_getreads(reads_file, genome_file, target, width)

def sim_genome_getreads(reads_file, genome_file, target, width):
  out_reads_file = 'extracted_sim_reads_c' + str(target) + '_s' + str(width) + '.fasta'
  out_genome_file = 'extracted_sim_genome_c' + str(target) + '_s' + str(width) + '.fasta'

  target_high = target + width / 2
  target_low = target - width / 2

  reads = []

  start = 0
  end = 0
  getdna = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if getdna:
        getdna = False
        # print line[target_low - start : end - target_high]
        reads.append(line[target_low - start : target_high - start])
      if line[0] == '>':
        start = line.split('/')[1].split('_')[1]
        end = line.split('/')[1].split('_')[2]
        start = int(start)
        end = int(end)
        if start < target_low and target_high < end:
          getdna = True
          reads.append(line.strip())

  with open(out_reads_file, 'w') as f:
    for line in reads:
      f.write(line + '\n')

  genome_section = []

  with open(genome_file) as f:
    for i, line in enumerate(f):
      if i == 0:
        genome_section.append(line.strip())
      if i == 1:
        genome_section.append(line[target_low:target_high])

  with open(out_genome_file, 'w') as f:
    for line in genome_section:
      f.write(line + '\n')

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start