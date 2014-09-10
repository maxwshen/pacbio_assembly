# Given a blasr sam output and its set of reads, orders the reads by genomic position

import sys
import string
import datetime
import random
import copy
import os
import read_fasta

def main():
  reads_file = sys.argv[1]  
  sam_file = sys.argv[2]
  order_reads(reads_file, sam_file)
  return

def order_reads(reads_file, sam_file):
  order = []    # List of tuples (header, pos)
  with open(sam_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      order.append(('/'.join(words[0].split('/')[:-1]), int(words[6])))

  order.sort(key = lambda o : o[1])

  new_reads = []
  h, r = read_fasta.read_fasta(reads_file)
  new_h = []
  for header in h:
    new_h.append(header[1:])
  for o in order:
    curr = o[0]
    print curr
    print new_h[0]
    new_reads.append(curr)
    new_reads.append(r[new_h.index(curr)])

  for line in new_reads:
    print line

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start