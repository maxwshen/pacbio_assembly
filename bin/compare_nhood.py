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

def main():
  nhood_file = sys.argv[1]
  reads_file = sys.argv[2]
  _k = int(sys.argv[3])
  _t = int(sys.argv[4])

  genome_file = '/home/mshen/research/data/e_coli_genome.fasta'
  
def compare_nhood(nhood_file):
  pass

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start