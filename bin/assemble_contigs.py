import sys
import string
import datetime
import random
import copy
import assembly
import locAL
import os

from collections import defaultdict

def main():
  # contigs_file format:
  # [start] [end] ACTG...
  #
  # highdegnodes_file format:
  # ACTG... Deg = X t = Y [in]correct
  contigs_file = sys.argv[1]
  highdegnodes_file = sys.argv[2] 
  genome_file = sys.argv[3]

  for i in range(12):
    _i = str(i)
    contigs_file = '/home/mshen/research/highcov__fold_s' + _i + '.t1.15.L0/contigsout.s' + _i + '.t1.15.L0.txt'
    highdegnodes_file = '/home/mshen/research/highcov__fold_s' + _i + '.t1.15.L0/highdegnodes.15.s' + _i + '.t1.kmers.out'
    genome_file = '/home/mshen/research/data/high_cov/ec_genome_rh_hc_n' + _i + '.fasta'

    supercontig = assemble_contigs(contigs_file, highdegnodes_file)
    checkAccuracy(supercontig, genome_file)

def assemble_contigs(contigs_file, highdegnodes_file):
  length_cutoff = 20
  contigs = [] # List of tuples (start, end, contig)
  with open(contigs_file) as f:
    for i, line in enumerate(f):
      if line[0] == '[' and line[:2] != '[]':
        line = line.translate(None, '[]')
        words = line.split()
        if len(words) == 3 and len(words[2]) > length_cutoff:
          contigs.append((int(words[0]), int(words[1]), words[2]))

  nodes = dict()  # Keys = kmers, Values = degree
  with open(highdegnodes_file) as f:
    for i, line in enumerate(f):
      kmer = line.split()[0]
      degree = int(line.split()[3])
      nodes[kmer] = degree
      lenkmer = len(kmer)

  while True:
    if find_overlap(contigs) != None:
      overlap, c1, c2 = find_overlap(contigs)
      c1kmer = c1[2]
      c2kmer = c2[2]
    else:
      break

    if overlap > lenkmer:
      print 'Overlap > lenkmer:', overlap, c1kmer, c2kmer
      sys.exit(0)

    print overlap
    print c1kmer[-lenkmer:], nodes[c1kmer[-lenkmer:]], c2kmer[:lenkmer], nodes[c2kmer[:lenkmer]]
    if nodes[c1kmer[-lenkmer:]] > nodes[c2kmer[:lenkmer]]:
      new_contig = c1kmer + c2kmer[overlap:]
    else:
      new_contig = c1kmer[:-overlap] + c2kmer

    print c1[0], c2[1], new_contig
    contigs.append((c1[0], c2[1], new_contig))
    del contigs[contigs.index(c1)]
    del contigs[contigs.index(c2)]

  return contigs[0][2]

def checkAccuracy(seq, genome_file):
  # Check the accuracy of the super contig
  with open(genome_file) as f:
    lines = f.readlines()
    genome = lines[1]

  (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(seq, genome, 1, -1, -1, -0.5)
  if alignLen != 0:
    score = float(matches) / float(alignLen)
  else:
    score = 0
    print 'ERROR: Divide by Zero! Alignment Length = 0'
  print 'Range:', bestxy[1] - alignLen, '-', bestxy[1]
  print 'Accuracy:', score
  print 'Error Rate:', 1 - score

  return score, alignLen

def find_overlap(contigs):
  if len(contigs) == 1:
    return None
  for i in range(len(contigs)):
    c1 = contigs[i]
    for c2 in [c for c in contigs if c != c1]:
      c1start, c1end, c1kmer = c1
      c2start, c2end, c2kmer = c2
      if c1start < c2start < c1end:
        # +1 because the ending position is 1 less than it should be. BUG to fix
        # +1 because there are 2 nucleotides between positions 1 and 2
        return c1end - c2start + 2, c1, c2
  return None

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start