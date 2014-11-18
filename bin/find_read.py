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
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  read = find_read('>' + header, reads_file)
  print read

def find_read(header, reads_file):
  get_line = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if get_line:
        return header + '\n' + line
      if line.strip() == header:
        get_line = True
  return None

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start