import sys
import copy
import string
import datetime
import random
import numpy as np

import locAL

from collections import defaultdict

def main():
  if len(sys.argv) != 3:
    print 'Usage: python pctAccuracy.py <reads.fasta> <genome.fasta>'
    sys.exit(0)

  pctAccuracy(sys.argv[1], sys.argv[2])


# bestseq1 will force all of seq1 to be in the alignment
def pctAccuracy(reads, genome):
  h, r = read_fasta(genome)
  genomeSeq = r[0]

  h, r = read_fasta(reads)
  accuracies = []
  for seq in r:
    (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(seq, genomeSeq, 1, -1, -1, -0.5)
    if alignLen > 0:
      accuracy = float(matches)/float(alignLen)
    else:
      accuracy = 0
      print 'WARNING: Alignment Length = 0'
    accuracies.append(accuracy)
    print 'Accuracy:', accuracy, '\tError Rate:', 1 - accuracy
    
  print accuracies
  print np.mean(accuracies)


def read_fasta(fasta_file):
  # Builds header and reads lists from a fasta file
  headers = []
  reads = []
  get_read = False
  with open(fasta_file) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] == '>' or line[0] == '@':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        headers.append(line.strip())
        get_read = True
      elif get_read:
        curr_read += line.strip()
  if curr_read != '':
    reads.append(curr_read)
  if len(headers) == 0 or len(reads) == 0:
    print 'ERROR: Bad fasta file', fasta_file
    sys.exit(0)
  return headers, reads

# main
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start
  main()
  end = datetime.datetime.now()
  print '\nEnd:', end, '\nTotal:', end - start
