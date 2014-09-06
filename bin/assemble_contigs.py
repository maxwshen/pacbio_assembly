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

  supercontigs = assemble_contigs(contigs_file, highdegnodes_file)
  for supercontig in supercontigs:
    checkAccuracy(supercontig[2], genome_file)
  return

  # Batch
  for i in range(12):
    _i = str(i)
    contigs_file = '/home/mshen/research/highdegassembly_stats/highcov__fold_s' + _i + '.t1.15.L0/contigsout.s' + _i + '.t1.15.L0.txt'
    highdegnodes_file = '/home/mshen/research/highdegassembly_stats/highcov__fold_s' + _i + '.t1.15.L0/highdegnodes.15.s' + _i + '.t1.kmers.out'
    genome_file = '/home/mshen/research/data/high_cov/ec_genome_rh_hc_n' + _i + '.fasta'

    supercontigs = assemble_contigs(contigs_file, highdegnodes_file)
    for supercontig in supercontigs:
      checkAccuracy(supercontig[2], genome_file)

def assemble_contigs(contigs_file, highdegnodes_file):
  length_cutoff = 18
  contigs = [] # List of tuples (start, end, contig)
  with open(contigs_file) as f:
    for i, line in enumerate(f):
      if line[0] == '[' and line[:2] != '[]':
        line = line.translate(None, '[]')
        words = line.split()
        if len(words) == 3 and len(words[2]) > length_cutoff:
          contigs.append((float(words[0]), float(words[1]), words[2]))

  nodes = dict()  # Keys = kmers, Values = degree
  with open(highdegnodes_file) as f:
    for i, line in enumerate(f):
      kmer = line.split()[0]
      degree = int(line.split()[3])
      _t = int(line.split()[6])
      nodes[kmer] = (degree, _t)
      lenkmer = len(kmer)

  while True:
    ans, overlap, c1, c2 = find_overlap(contigs, lenkmer)
    if ans != False:
      c1kmer = c1[2]
      c2kmer = c2[2]
    else:
      break

    if overlap > lenkmer:
      print 'Overlap > lenkmer:', overlap, c1, c2
      sys.exit(0)

    print overlap
    # print c1kmer[-lenkmer:], nodes[c1kmer[-lenkmer:]], c2kmer[:lenkmer], nodes[c2kmer[:lenkmer]]
    # First by t, then by degree if t's are equal
    if c1kmer[-lenkmer:] not in nodes.keys():
      print c1kmer[-lenkmer:], 'not in high-deg nodes'
      c1deg = 0
      c1t = 0
    else:
      c1deg = nodes[c1kmer[-lenkmer:]][0]
      c1t = nodes[c1kmer[-lenkmer:]][1]
    if c2kmer[:lenkmer] not in nodes:
      print c2kmer[:lenkmer], 'not in high-deg nodes'
      c2deg = 0
      c2t = 0
    else:
      c2deg = nodes[c2kmer[:lenkmer]][0]
      c2t = nodes[c2kmer[:lenkmer]][1]
      
    if c1t > c2t:
      new_contig = c1kmer + c2kmer[overlap:]
    elif c1t < c2t:
      new_contig = c1kmer[:-overlap] + c2kmer
    else:
      if c1deg > c2deg:
        new_contig = c1kmer + c2kmer[overlap:]
      else:
        new_contig = c1kmer[:-overlap] + c2kmer

    print c1, c2
    print c1[2] + '-' * (len(c2[2]) - overlap)
    print '-' * (len(c1[2]) - overlap) + c2[2]
    print new_contig 
    contigs.append((c1[0], c2[1], new_contig))
    del contigs[contigs.index(c1)]
    del contigs[contigs.index(c2)]

  return contigs

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

def find_overlap(contigs, lenkmer):
  if len(contigs) == 1:
    return (False, -1, None, None)

  # Remove supercontigs inside of others
  remove_set = set()
  for i in range(len(contigs)):
    c1 = contigs[i]
    for c2 in [c for c in contigs if c != c1]:
      c1start, c1end, c1kmer = c1
      c2start, c2end, c2kmer = c2
      if c1start < c2start < c2end < c1end:
        remove_set.add(c2)
  for r in remove_set:
    del contigs[contigs.index(r)]   

  for i in range(len(contigs)):
    c1 = contigs[i]
    for c2 in [c for c in contigs if c != c1]:
      c1start, c1end, c1kmer = c1
      c2start, c2end, c2kmer = c2
      if c1start < c2start < c1end < c2end:
        # +1 because the ending position is 1 less than it should be. BUG to fix
        # +1 because there are 2 nucleotides between positions 1 and 2
        if int(c1end - c2start + 1) <= lenkmer:
          return True, int(c1end - c2start + 1), c1, c2   
  return (False, -1, None, None)

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start