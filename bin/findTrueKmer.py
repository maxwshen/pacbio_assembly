# findTrueKmer.py
# 
# Takes in reads. Attempts to find "true" kmers by 
# creating a graph where nodes are kmers and edges connect
# kmers that are 1 insertion or 1 deletion away from each
# other.
#
# Implementated with a table approach
#
# Genome file is used to measure accuracy.

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import GenerateIndelKmers
import locAL

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])
  cutoff = int(sys.argv[4])
  t_atleast = int(sys.argv[5])
  if sys.argv[6] == 'True':
    fn = True
  else:
    fn = False

  findTrueKmer(reads_file, genome_file, _k, cutoff, t_atleast, fn)
  return 

  for t_atleast in range(1, 11):
    findTrueKmer(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), t_atleast)


def generateRandomKmer(_k):
  nt = ['A', 'C', 'T', 'G']
  kmer = ''
  for i in range(_k):
    r = int(np.random.randint(4, size=1))
    kmer += nt[r]
  return kmer

def graphviz(genomeKmers, degrees, kmers, cutoff):
  # genomeKmers is a set of all kmers in the reads
  # degrees is a dict. Keys = kmer, Values = degree
  # kmers is a dict. Keys = kmer, Values = t

  numToOutput = cutoff
  num = copy.copy(numToOutput)
  with open('findTrueKmers.gv', 'w') as f:
    f.write('digraph G {\n')
    for key in sorted(degrees, key=degrees.get, reverse=True):
      if num == 0:
        break
      if key in genomeKmers:
        f.write(key[:-1] + ' -> ' + key[1:] + ';\n')
      else:
        f.write(key[:-1] + ' -> ' + key[1:] + ' [color="red"];\n')
      num -= 1
    f.write('}')


def filter_neighbors(degrees, cutoff):
  # Input: Dict, keys = kmers, values = degrees
  new_degrees = dict()
  removed = set()
  count = 0
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if key not in removed:
      new_degrees[key] = degrees[key]
      neighbors = find_neighbors(degrees, key)
      for n in [i for i in neighbors if i in degrees]:
        removed.add(n)
    else:
      count += 1
    if len(new_degrees) > cutoff:
      break

  print 'Filtered', count, 'kmers'
  return new_degrees

def find_neighbors(degrees, kmer):
  del_kmers = GenerateIndelKmers.genDelKmers(kmer, 1)[1]
  remove_list = []
  for dk in del_kmers:
    if dk == kmer[:-1] or dk == kmer[1:]:
      remove_list.append(dk)
  for dk in remove_list:
    del_kmers.remove(dk)
  remove_list = []
  ins_kmers = GenerateIndelKmers.genInsKmers(kmer, 1)[1]
  for ik in ins_kmers:
    if kmer == ik[:-1] or kmer == ik[1:]:
      remove_list.append(ik)
  for ik in remove_list:
    ins_kmers.remove(ik)

  del_ins = set()
  for dk in del_kmers:
    for i in GenerateIndelKmers.genInsKmers(dk, 1)[1]:
      del_ins.add(i)
  for ik in ins_kmers:
    for i in GenerateIndelKmers.genDelKmers(ik, 1)[1]:
      del_ins.add(i)

  del_ins.remove(kmer)

  neighbors = [i for i in del_ins if i in degrees]
  return neighbors

def consensus(reads):
  # Constructs a multiple alignment where each column has exactly 1 kind of nt
  allnt = [list(reads[0])]
  mat = [list(reads[0])]
  for r in reads[1:]:
    print mat
    curr = []
    for i in range(len(r)):
      char = r[i]
      if char != allnt[i]:
        allnt = allnt[0:i] + [char] + allnt[i:]
        for j in range(len(mat)):
          mat[j] = mat[j][0:i] + ['-'] + mat[j][i:]  
      curr.append(char)
    mat.append(curr)

  for m in mat:
    print m


def findTrueKmer(reads_file, genome_file, _k, cutoff, t_cutoff, fn):
  _kplus = _k + 1
  _kminus = _k - 1
  isdna = False
  readcount = 0
  kmers = dict()
  kplusmers = dict()
  kminusmers = dict()

  with open(reads_file) as f:
    for i, line in enumerate(f):
      if isdna:
        isdna = False
        dna = line.strip()
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j+_k]
          if kmer in kmers:
            kmers[kmer] += 1
          else:
            kmers[kmer] = 1
        for j in range(len(dna) - _kminus + 1):
          kmer = dna[j:j+_kminus]
          if kmer in kminusmers:
            kminusmers[kmer] += 1
          else:
            kminusmers[kmer] = 1
        for j in range(len(dna) - _kplus + 1):
          kmer = dna[j:j+_kplus]
          if kmer in kplusmers:
            kplusmers[kmer] += 1
          else:
            kplusmers[kmer] = 1
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        isdna = True

  # Generate kmers from the genome
  with open(genome_file) as f:
    lines = f.readlines()
    genome = lines[1].strip()
    middleGenome = genome[100:900]
  genomeKmers = set()
  genomeMidKmers = set()
  for i in range(len(middleGenome) - _k + 1):
    genomeMidKmers.add(middleGenome[i:i+_k])
  for i in range(len(genome) - _k + 1):
    genomeKmers.add(genome[i:i+_k])

  for kmer in genomeKmers:
    if kmer not in kmers:
      kmers[kmer] = 0

  # Find degrees of all kmers
  degrees = dict()    # Key = kmer, Value = degree
  for kmer in kmers:
    degree = 0
    for del_kmer in GenerateIndelKmers.genDelKmers(kmer, 1)[1]:
      if del_kmer in kminusmers:
        if del_kmer != kmer[:-1] and del_kmer != kmer[1:]:
          # degree += 1
          degree += kminusmers[del_kmer]

    for ins_kmer in GenerateIndelKmers.genInsKmers(kmer, 1)[1]:
      if ins_kmer in kplusmers:
        if kmer != ins_kmer[:-1] and kmer != ins_kmer[1:]:
          degree += kplusmers[ins_kmer]

    degrees[kmer] = degree + kmers[kmer]

  # Generate degrees by finding all possible deletions of (k+1)-mers
  # for kplusmer in kplusmers:
  #   normal_kmers = GenerateIndelKmers.genDelKmers(kplusmer, 1)[1]
  #   for kmer in normal_kmers:
  #     if kmer != kplusmer[:-1] and kmer != kplusmer[1:]:
  #       if kmer in degrees:
  #         degrees[kmer] += kplusmers[kplusmer]
  #       else:
  #         degrees[kmer] = kplusmers[kplusmer]
  #       if kmer not in kmers:
  #         kmers[kmer] = 0

  for kmer in degrees:
    degrees[kmer] += kmers[kmer]
  for kmer in kmers:
    if kmer not in degrees:
      degrees[kmer] = kmers[kmer]

  for n in find_neighbors(degrees, 'ATGCACTGGGCATAC'):
    print n, degrees[n]

  consensus(find_neighbors(degrees, 'ATGCACTGGGCATAC'))
  return

  if fn:
    degrees = filter_neighbors(degrees, cutoff)

  numToOutput = cutoff
  numincorrect = 0
  current_deg = 0
  num = copy.copy(numToOutput)
  best = set()
  numoft = dict()
  for key in sorted(degrees, key=degrees.get, reverse=True):
    if kmers[key] >= t_cutoff:      # Filter by t
      if num >= 0:
        current_deg = degrees[key]
      elif degrees[key] != current_deg:
        break

      if key in genomeKmers:
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'correct'
      else:
        print key, 'Deg =', degrees[key], 't =', kmers[key], 'incorrect'
        numincorrect += 1
      if kmers[key] in numoft:
        numoft[kmers[key]] += 1
      else:
        numoft[kmers[key]] = 1
      best.add(key)
      num -= 1


  # print '\nCorrect Only:'
  # for key in sorted(degrees, key=degrees.get, reverse=True):
  #   if key in genomeKmers:
  #     print key, 'Deg =', degrees[key], 't =', kmers[key], 'correct'  
  #     genomeKmers.remove(key)
  # print 'Missing:', genomeKmers

  graphviz(genomeKmers, degrees, kmers, cutoff)

  # print '\nIncorrect'
  # for key in sorted(degrees, key=degrees.get, reverse=True):
  #   if key not in genomeKmers:
  #     print key, 'Deg =', degrees[key], 't =', kmers[key]

  print 't-cutoff:', t_cutoff, 'Incorrect:', numincorrect, 'Total:', cutoff - num
  return

  # print numoft
  total = float(numToOutput)
  t2orhigher = 0
  t3orhigher = 0
  t4orhigher = 0
  for key in numoft.keys():
    if key >= 2:
      t2orhigher += numoft[key]
  t3orhigher = t2orhigher - numoft[2]
  t4orhigher = t3orhigher - numoft[3]

  # Generate incorrect kmers that do not exist in genome
  numIncorrectKmers = 5000
  genomeKmers = set()
  for i in range(len(genome) - _k + 1):
    genomeKmers.add(genome[i:i+_k])
  wrongKmers = set()
  while len(wrongKmers) < numIncorrectKmers:
    randKmer = generateRandomKmer(_k)
    if randKmer not in genomeKmers:
      wrongKmers.add(randKmer)

  # Find degree of incorrect kmers
  degreesWrong = dict()
  for kmer in wrongKmers:
    del_kmers = GenerateIndelKmers.genDelKmers(kmer, 1)[1]
    degree = 0
    for del_kmer in del_kmers:
      if del_kmer in kminusmers:
        degree += 1

    ins_kmers = GenerateIndelKmers.genInsKmers(kmer, 1)[1]
    for ins_kmer in ins_kmers:
      if ins_kmer in kplusmers:
        degree += 1 

    degreesWrong[kmer] = degree

  print '\n', numIncorrectKmers, 'Incorrect kmers, sorted\n', 
  for key in sorted(degreesWrong, key=degrees.get, reverse=True):
    # print key, 'Deg =', degreesWrong[key]
    pass

  # return

  print _k

  numPerfect = 0
  matchScores = []
  perfectStarts = []
  perfectLens = []
  for kmer in best:
    (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(kmer, genome, 1, -1, -1, -0.5)
    # print 'Range:', bestxy[1] - alignLen, '-', bestxy[1]
    score = float(matches)/float(alignLen)
    matchScores.append(score)
    
    if score == 1:
      # Only consider correct kmers as within center 800bp
      if 100 < bestxy[1] - alignLen < 900: 
        print kmer, bestxy[1] - alignLen
        numPerfect += 1
        perfectStarts.append(bestxy[1] - alignLen)
        perfectLens.append(alignLen)
  perfectLens = [x for (y,x) in sorted(zip(perfectStarts, perfectLens))]
  perfectStarts = sorted(perfectStarts)

  # for score in matchScores:
  #   print str(score) + '\t',
  # print '\n'

  perfectRanges = [] # stores tuples
  start = -1
  end = -1
  for i in range(len(genome)):
    if i in perfectStarts:
      if end < i + perfectLens[perfectStarts.index(i)]:
        end = i + perfectLens[perfectStarts.index(i)]
      if start == -1:
        start = i
    if i == end:
      perfectRanges.append((start, end))
      start = -1
      end = -1

  totalCovered = 0
  for (x, y) in perfectRanges:
    totalCovered += y - x

  # print '2 or higher t:', t2orhigher, float(t2orhigher)/total
  # print '3 or higher t:', t3orhigher, float(t3orhigher)/total
  # print '4 or higher t:', t4orhigher, float(t4orhigher)/total
  # print '# kmers:', numToOutput
  # print '# perfect:', numPerfect, '\t', float(numPerfect*100)/float(numToOutput), '%'
  # print 'Avg. accuracy:', 100*sum(matchScores)/len(matchScores), '%'
  # print 'Covered ranges:', perfectRanges
  # print 'Total covered:', totalCovered, '\tPercent:', float(totalCovered)/float(len(genome))


if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start