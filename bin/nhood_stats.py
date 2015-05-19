# Finds the specificity and sensitivity of a neighborhood compared to the ground truth
# Ground truth = genome_nhood.txt, produced by blasr alignment of all reads to the genome
#   and taking any that aligns >7kb with the region where a read aligns.

import sys, string, datetime, random, copy, os, commands, fnmatch
from collections import defaultdict

ground_truth = 'genome_nhood.txt'
inp = 'nhoods_temp_creads.outrx_27_6_rc_v2.out'

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
    nhood[base] = rest
    if base not in nhood[base]:
      nhood[base].append(base)

specificity = []
sensitivity = []
for key in nhood.keys():
  if key not in true_nhoods:
    specificity.append(0)
  if key in true_nhoods:
    intersect = len(set(nhood[key]).intersection(true_nhoods[key]))
    sensitivity.append(float(intersect) / float(len(true_nhoods[key])))
    specificity.append(float(intersect) / float(len(nhood[key])))

print 'sensitivity:', float(sum(sensitivity)) / float(len(sensitivity))
print 'specificity:', float(sum(specificity)) / float(len(specificity))
