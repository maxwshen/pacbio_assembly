# spec_multi_align.py
# Max Shen - maxwshen@gmail.com
#
# Performs a special multiple alignment where each column
# has exactly one nucleotide in it.
# Used to find the consensus between neighbors in the
# Hammer graph.

import sys
import string
import datetime
import random
import copy
import os
import numpy as np

import read_fasta as rf

from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  h, r = rf.read_fasta(reads_file)

  spec_multi_align(r)
  return 


def spec_multi_align(seqs):
  # Input: List of DNA seqs, at least 2

  s1, s2 = align(seqs[0], seqs[1])
  alignment = [list(s1), list(s2)]
  master = get_master(s1, s2)

  for seq in seqs[2:]:
    mas, curr_seq = align(''.join(master), seq)
    alignment = adjust_alignment(alignment, mas)
    alignment.append(list(curr_seq))
    master = get_master(mas, curr_seq)

  final_align = [''.join(a) for a in alignment]
  consensus = ''
  for i in range(len(final_align[0])):
    dashes = 0
    nt_count = 0
    nt = ''
    for row in final_align:
      if row[i] == '-':
        dashes += 1
      else:
        nt = row[i]
        nt_count += 1
    if nt_count >= dashes:
      consensus += nt
  
  for a in final_align:
    print a
  print consensus, len(consensus)
  return consensus

def adjust_alignment(alignment, mas):
  for i in range(len(mas)):
    if mas[i] == '-':
      for row in alignment:
        row.insert(i, '-')
  return alignment

def get_master(seq1, seq2):
  master = []
  for i in range(len(seq1)):
    if seq1[i] != '-':
      master.append(seq1[i])
    else:
      master.append(seq2[i])
  return master

def align(q1, q2):  
  match = 1
  indel = -1
  q1 = ' ' + q1
  q2 = ' ' + q2
  s = [[0]*(len(q1)) for i in range(len(q2))]
  bt = [['']*len(q1) for i in range(len(q2))]

  # Base case
  for i in range(1, len(q2)):
    s[i][0] = s[i-1][0] + indel
    bt[i][0] = 'up'
  for j in range(1, len(q1)):
    s[0][j] = s[0][j-1] + indel
    bt[0][j] = 'left'
  

  for i in range(1, len(q2)):
    # print i
    for j in range(1, len(q1)):
      ups = s[i-1][j] + indel
      lefts = s[i][j-1] + indel
      if q2[i] == q1[j]:
        diags = s[i-1][j-1] + match
      else:
        diags = -100
      s[i][j] = max(ups, lefts, diags)

      # Update bt
      if s[i][j] == diags:
        bt[i][j] = 'diag'
      elif s[i][j] == lefts:
        bt[i][j] = 'left'
      elif s[i][j] == ups:
        bt[i][j] = 'up'

  best = {'i': len(q2) - 1, 'j': len(q1) - 1, 'score': s[len(q2) - 1][len(q1) - 1]}

  # print s
  # for i in bt:
    # print i
  f1 = []
  f2 = []
  seq1, seq2 = printBT(s, bt, q1, q2, best['i'], best['j'])
  # print 'Score =', s[best['i']][best['j']]

  return seq1, seq2

def printBT(s, bt, q1, q2, i, j):
  f1 = []
  f2 = []
  while True:
    # print s[i][j]
    if i == 0 and j == 0:
      break
    elif bt[i][j] == 'up':
      f1.insert(0, '-')
      f2.insert(0, q2[i])
      i -= 1
    elif bt[i][j] == 'left':
      f1.insert(0, q1[j])
      f2.insert(0, '-')
      j -= 1
    else:
      f1.insert(0, q1[j])
      f2.insert(0, q2[i])
      i -= 1
      j -= 1

  seq1 = ''
  mid = ''
  seq2 = ''
  matches = float(0)
  posmatches = float(0)
  mismatches = 0
  gapopen = 0
  gapextend = 0
  totalLen = float(len(f1))
  for i in range(len(f1)):
    seq1 += f1[i]
    seq2 += f2[i]
    if f1[i] == f2[i]:
      mid += '|'
      matches += 1
    elif f1[i] == '-' or f2[i] == '-':
      mid += ' '
      if i > 0:
        if f1[i-1] == '-' or f2[i-1] == '-':
          gapextend += 1
        else:
          gapopen += 1

  # print seq1, '\n', mid, '\n', seq2
  # print 'Accuracy =', 100*(matches+posmatches) / totalLen, '%'
  # print 'Length =', len(f1)
  # print 'M:', int(matches), '\tPM:', int(posmatches), '\tMM:', mismatches, '\tGO:', gapopen, '\tGE:', gapextend
  
  return seq1, seq2

if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start