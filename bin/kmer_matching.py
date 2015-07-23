# Given a single error-corrected read that hopefully corresponds to some genomic region,
# try to find reads that match this genomic region based on shared k-mers

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf

from collections import defaultdict

def main():
  ec_seq_file = sys.argv[1]
  _k = int(sys.argv[2])
  cutoff = int(sys.argv[3])
  read_file = sys.argv[4]
  # read_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  out_file = 'km_' + ec_seq_file.translate(None, '/') + '.fasta'
  kmer_matching(ec_seq_file, read_file, _k, cutoff, file_bool = True)
  return

  # Batch
  # fold = '/home/mshen/research/yu/N250_4/'
  fold = '/home/mshen/research/yu_ec_22.4_500_nhoods/'
  for name in os.listdir(fold):
    ec_seq_file = name
    out_file = 'km_' + ec_seq_file.translate(None, '/') + '.fasta'    
    kmer_matching_and_checking(fold + ec_seq_file, read_file, _k, cutoff, out_file)
  return

def kmer_matching(ec_seq_file, read_file, _k, cutoff, file_bool = True):
  # ec_seq_file should be a fasta file with only one sequence
  # if file is false, ec_seq_file is the sequence directly.

  if file_bool:
    hs, rs = rf.read_fasta(ec_seq_file)
    ec_seq = rs[0]
  else:
    ec_seq = ec_seq_file

  print 'Finding k-mers in base read...'
  kmers = set()
  for i in range(len(ec_seq) - _k + 1):
    kmers.add(ec_seq[i:i + _k])
  print 'Found', len(kmers), 'kmers'

  print 'Finding k-mers in all other reads...'
  reads = dict()    # Key = header, value = num shared kmers
  hr, rr = rf.read_fasta(read_file)
  for i in range(len(rr)):
    if i % 1000 == 0:
      print i
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  filtered_heads = []
  for key in sorted(reads, key = reads.get, reverse = True):
    # print key, reads[key]
    if reads[key] < cutoff or len(filtered_heads) > 50:
      break
    filtered_heads.append(key)

  return filtered_heads

def kmer_matching_and_checking(ec_seq_file, read_file, _k, cutoff, out_file):
  headers = kmer_matching(ec_seq_file, read_file, _k, cutoff)

  print '\n' + ec_seq_file + '\n' + str(len(headers))
  get_reads_from_headers(headers, read_file, out_file)
  return
  

  # Use blasr to check accuracy of kmer-matching found reads

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