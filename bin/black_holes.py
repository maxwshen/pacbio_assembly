# SAMPLE INPUT
# /home/mshen/research/e_coli_nh_ec_22.4_unwhole_100k/nhood_nh1077087_CAGTTTGCCCTTCCTGTTGCGCC.fasta_Cov5.fasta
# 14
# >m120114_011938_42177_c100247042550000001523002504251220_s1_p0/50279/642_8506/0_7864
# 166 1076320 1083736 7416

# Given a single error-corrected read that hopefully corresponds to some genomic region,
# try to find reads that match this genomic region based on shared k-mers

import sys
import string
import datetime
import random
import copy
import os
import commands
import numpy as np

from collections import defaultdict

def main():
  km_out_file = sys.argv[1]

  black_holes(km_out_file)
  return

def black_holes(km_out_file):
  grab_pos = False
  currpos = None
  pos = [0 for i in range(100001)]
  collected_headers = set()

  total = 0
  with open(km_out_file) as f:
    for i, line in enumerate(f):
      if grab_pos:
        grab_pos = False
        small = int(line.split()[1])
        large = int(line.split()[2])
        if small < currpos < large:
          # print large, small, large - small
          total += large - small
          for j in range(small, min(large + 1, 1100000)):
            pos[j - 1000000] += 1
      if line[0] == '/':
        currpos = int(line.split('_')[7][2:])
      if line[0] == '>':
        if line not in collected_headers:
          collected_headers.add(line)
          grab_pos = True


  above = [0 for i in range(30)]
  for i in pos:
    for j in range(len(above)):
      if i >= j:
        above[j] += 1

  print 'input:', km_out_file
  print above
  print 'min cov', min(pos), ', max cov', max(pos)
  print len(collected_headers), 'unique reads found'
  print total, 'total bp found'
  print float(total) / float(len(collected_headers)), 'avg read len'
  print float(total) / float(100000), 'avg cov'
  print pos.index(min(pos)), 'is a position of lowest cov'


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start