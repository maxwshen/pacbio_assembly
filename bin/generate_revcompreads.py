import sys, string, datetime, random, copy, os
import numpy as np
from collections import defaultdict
import find_read

def main():
  # reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  for num in range(5, 60, 5):
    for samp in range(1, 4):
      snum = str(num)
      ssamp = str(samp)
      reads_file = '/home/jeyuan/pacbio_undersample/new_reads/reads.20k.cov' + snum + '.samp' + ssamp + '.fasta'
      new_reads_file = '/home/mshen/research/data/undersampled_20k_rc/reads.20k.' + snum + 'x.' + ssamp + 'rc.fasta'
      generate_revcompreads(reads_file, new_reads_file)
  return


def generate_revcompreads(reads_file, new_reads_file):
  mapper = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
  rev_prefix = '>R_'
  with open(reads_file) as f:
    with open(new_reads_file, 'w') as fout:
      dna = ''
      prev_header = ''
      for i, line in enumerate(f):
        if i % 100000 == 0:
          print i
        if i == 0:
          prev_header = line.strip()
        if i > 0:
          if line[0] == '>':
            fout.write(rev_prefix + prev_header[1:] + '\n')
            fout.write(reverse_complement(dna, mapper) + '\n')
            fout.write(prev_header + '\n')
            fout.write(dna + '\n')
            dna = ''
            prev_header = line.strip()
          else:
            dna += line.strip()
      fout.write(rev_prefix + prev_header[1:])
      fout.write(reverse_complement(dna, mapper))
      fout.write(prev_header)
      fout.write(dna)
  return    

def reverse_complement(dna, mapper):
  new_dna = ''
  for nt in dna[::-1]:
    new_dna += mapper[nt]
  return new_dna

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start