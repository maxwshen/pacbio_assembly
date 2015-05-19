# Finds the specificity and sensitivity of a neighborhood compared to the ground truth
# Ground truth = genome_nhood.txt, produced by blasr alignment of all reads to the genome
#   and taking any that aligns >7kb with the region where a read aligns.

import sys, string, datetime, random, copy, os, commands, fnmatch
from collections import defaultdict

ground_truth = '/home/mshen/research/data/genome_nhood.txt'
# ground_truth = '/home/mshen/research/data/genome_nhood_all.txt'
inp = '/home/mshen/research/data/nhoods_temp_creads.outrx_27_6_rc_v2.out'

true_nhoods = defaultdict(list)
with open(ground_truth) as f:
  for i, line in enumerate(f):
    base = line.split(':')[0]
    rest = ' '.join(line.split(':')[1:]).split()
    true_nhoods[base] = rest

nhood = defaultdict(list)
with open(inp) as f:
  for i, line in enumerate(f):
    base = line.split()[0]
    rest = line.split()[1:]
    if len(rest) > 0:
      nhood[base] = rest
      if base not in nhood[base]:
        nhood[base].append(base)

specificity = []
sensitivity = []
num_false_nhood = 0
for key in nhood.keys():
  if key not in true_nhoods:
    # specificity.append(0)
    # print key, 'has no ground truth nhood'
    num_false_nhood += 1
    # print key, nhood[key]
  if key in true_nhoods:
    intersect = len(set(nhood[key]).intersection(true_nhoods[key]))
    sensitivity.append(float(intersect) / float(len(true_nhoods[key])))
    specificity.append(float(intersect) / float(len(nhood[key])))

print 'sensitivity:', float(sum(sensitivity)) / float(len(sensitivity))
print 'specificity:', float(sum(specificity)) / float(len(specificity))
print 'found', num_false_nhood, 'number of false nhoods (where no ground truth nhood exists)'

# print 'sensitivity list:', '\n'.join(str(s) for s in sensitivity)
# print 'specificity list:', '\n'.join(str(s) for s in specificity)
