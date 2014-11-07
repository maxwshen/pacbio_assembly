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

import read_fasta as rf

from collections import defaultdict

def main():
  ec_seq_file = sys.argv[1]
  _k = int(sys.argv[2])
  cutoff = int(sys.argv[3])
  read_file = sys.argv[4]
  # read_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'

  kmer_matching(ec_seq_file, read_file, _k, cutoff)
  return

  # Batch
  fold = '/home/mshen/research/yu/N250_4/'
  for name in os.listdir(fold):
    ec_seq_file = name
    kmer_matching(fold + ec_seq_file, read_file, _k, cutoff)
  return

def kmer_matching(ec_seq_file, read_file, _k, cutoff):
  # ec_seq_file should be a fasta file with only one sequence

  hs, rs = rf.read_fasta(ec_seq_file)
  ec_seq = rs[0]

  kmers = set()
  for i in range(len(ec_seq) - _k + 1):
    kmers.add(ec_seq[i:i + _k])

  reads = dict()    # Key = header, value = num shared kmers
  hr, rr = rf.read_fasta(read_file)
  for i in range(len(rr)):
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  headers = []
  for key in sorted(reads, key = reads.get, reverse = True):
    if reads[key] < cutoff:
      break
    headers.append(key)

  print '\n' + ec_seq_file + '\n' + str(len(headers))

  new_ec_seq = ec_seq_file.translate(None, '/')
  out_file = 'km_' + new_ec_seq + '.fasta'
  get_reads_from_headers(headers, read_file, out_file)

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  blasr_options = '-bestn 1'
  blasr_out = commands.getstatusoutput(blasr_exe + ' ' + out_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]

  to_print = []
  for line in blasr_out.splitlines():
    head = '>' + '/'.join(line.split()[0].split('/')[:-1])
    start_pos = line.split()[6]
    end_pos = line.split()[7]
    length = int(end_pos) - int(start_pos)
    info = (head, str(reads[head]), start_pos, end_pos, str(length))
    to_print.append(info)

  for t in sorted(to_print, key = lambda tup: int(tup[1]), reverse = True):
    print t[0]
    print '\t'.join(t[1:])

  return

def get_reads_from_headers(headers, read_file, out_file):
  # headers_file should have exact full headers as the first "word" on each line 
  reads = ''
  h, r = rf.read_fasta(read_file)
  for i in range(len(h)):
    if h[i] in headers:
      reads += h[i] + '\n' + r[i] + '\n'
  with open(out_file, 'w') as f:
    f.write(reads)
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start