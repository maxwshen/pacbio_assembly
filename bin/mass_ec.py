import sys, string, datetime, random, copy, os, commands, fnmatch
import numpy as np

import read_fasta as rf
import find_read
from collections import defaultdict

def main():
  nhood_fold = '/home/mshen/research/e_coli_nhoods_500_22.4_unwhole'
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  mass_ec(nhood_fold, reads_file)

def mass_ec(fold, reads_file):
  ec_tool = '/home/mshen/research/bin/read_correction.sh'

  for name in os.listdir(fold):
    headers = []
    with open(fold + '/' + name) as f:
      for i, line in enumerate(f):
        if line[0] == '>':
          headers.append(line.strip())

    fullreads = find_read.find_reads(headers, reads_file)
    full_reads_file = name + 'temp_full_reads.fasta'
    with open(full_reads_file, 'w') as f:
      f.write(fullreads)

    status = commands.getstatusoutput(ec_tool + ' ' + full_reads_file)[1]

    considered = []
    for fn in os.listdir('.'):
      if fnmatch.fnmatch(fn, name + '*_read_*'):
        considered.append(fn)

    best_file = ''
    best_len = 0
    for fn in considered:
      with open(fn) as f:
        read = f.readlines()[1]
        if len(read) > best_len:
          best_len = len(read)
          best_file = fn

    status = commands.getstatusoutput('mv ' + best_file + ' ' + 'best_' + name + '.fasta')[1]
    status = commands.getstatusoutput('rm -rf *_read_*')[1]
    status = commands.getstatusoutput('rm -rf *_full_*')[1]



if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start