# Performs iterative error correction of reads.
# Beginning with an arbitrary error corrected consensus,
#   1. Use k-mer matching to find corresponding reads
#   2. Use blasr to align the EC consensus to these reads, finding indices for replacement
#   3. Replace EC consensus into reads, then grab next 500bp and perform error correction
# This process is repeated until the coverage falls too low for proper error correction

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf
import find_read
import kmer_matching
from collections import defaultdict

def main():
  reads_file = '/home/mshen/research/data/e_coli_genome.fasta'
  ec_reads_file = '/home/mshen/research/data/EC_e_coli_genome.fasta'
  commands.getstatusoutput('cp ' + reads_file + ' ' + new_reads_file)
  consensus_file = sys.argv[1]

  iterative_ec(ec_reads_file, ec_tool)

def iterative_ec(consensus_file, ec_reads_file):
  ec_tool = '/home/mshen/research/bin/error_correction.sh'

  _k = 15
  cutoff = 5
  kmer_matching.kmer_matching(consensus_file, ec_reads_file, 15, 5)

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start