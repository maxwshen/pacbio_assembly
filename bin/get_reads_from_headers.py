# Given a single error-corrected read that hopefully corresponds to some genomic region,
# try to find reads that match this genomic region based on shared k-mers

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import read_fasta as rf

from collections import defaultdict

def main():
  headers_file = sys.argv[1]
  read_file = sys.argv[2]

  get_reads_from_headers(headers_file, read_file)
  return

def get_reads_from_headers(headers_file, read_file):
  # headers_file should have exact full headers as the first "word" on each line 
  headers = set()
  with open(headers_file) as f:
    for i, line in enumerate(f):
      headers.add(line.split()[0].strip())

  reads = ''
  h, r = rf.read_fasta(read_file)
  for i in range(len(h)):
    if h[i] in headers:
      reads += h[i] + '\n' + r[i] + '\n'

  print reads
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start