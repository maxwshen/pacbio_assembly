import sys
import string
import datetime
import random
import copy
import os
import commands

from collections import defaultdict

def main():
  header = sys.argv[1]
  # reads_file = sys.argv[2]
  # reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  # reads_file = '/home/mshen/research/NEWREADS_22.4_rmhomo.fasta'
  # reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  reads_file = '/home/mshen/research/data/reads.20k.rc.fasta'
  read = find_read('>' + header, reads_file)
  print read

def find_read(header, reads_file):
  get_line = False
  read = ''
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if get_line:
        if line[0] == '>':
          break
        read += line.strip()
      if line[0] == '>':
        curr_header = line.strip()
        curr_header = ''.join(line.strip().split()[:-1])
        if curr_header == header:
          get_line = True
  if len(read) > 0:
    return header + '\n' + read
  return None

def find_reads(headers, reads_file):
  get_line = False
  output = ''
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if get_line:
        output += header + '\n' + line
        get_line = False
      if line.strip() in headers:
        header = line.strip()
        get_line = True
  return output

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start