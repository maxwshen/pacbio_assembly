# locAL.py
# BIMM 185
# Max Shen A10082759
#
# python locAL.py <seq file1> <seq file2> -m <match> -s <mismatch>
# -go <gap-open> -ge gap-extend
#
# locAL gives the best local alignment of two seq files based on an affine
# cost function provided by the user, using the Smith-Waterman Algorithm.
# It outputs the best score, the resulting sequence, and its length.
#
# locAL is not space efficient.
#
# You can call locAL from another python script by 'import locAL'
# and using external(seq1, seq2, mscore, mmscore, goscore, gescore)
# It will return a tuple (alignLen, matches, mismatches, numgaps, numGapExtends)
#
# locAL also supports "fitting" alignment, that is, forcing all of seq1
# to be included in the alignment. 
# external_bestseq1(seq1, seq2, mscore, mmscore, goscore, gescore)
#
# Suggested Use:
# (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(seq, genomeSeq, 1, -1, -1, -0.5)
#  accuracy = float(matches)/float(alignLen)
#

import sys
import copy
import string
import datetime
import random

from collections import defaultdict

# Backtracking variable names
bt_zero = 0
bt_up = 1
bt_left = 2
bt_diag = 3

def main():
  if len(sys.argv) < 3 or len(sys.argv) % 2 == 0:
    print 'Usage: python locAL.py <seq file1> <seq file2>'
    print 'Optional Flags: -m <match> -s <mismatch> -go <gap-open> -ge <gap-extend>'
    sys.exit(0)

  # Default values
  ms = 1
  mms = -1
  go = -1
  ge = -0.5

  h1, r1 = read_fasta(sys.argv[1])
  h2, r2 = read_fasta(sys.argv[2])
  seq1 = r1[0]
  seq2 = r2[0]

  # Detect flags and change parameters
  if len(sys.argv) > 3:
    for param in sys.argv[3:]:
      if param == '-m':
        ms = float(sys.argv[sys.argv.index('-m') + 1])
      if param == '-s':
        mms = float(sys.argv[sys.argv.index('-s') + 1])
      if param == '-go':
        go = float(sys.argv[sys.argv.index('-go') + 1])
      if param == '-ge':
        ge = float(sys.argv[sys.argv.index('-ge') + 1])

  print 'Match = ', ms, '\tMismatch = ', mms
  print 'Gap Open = ', go, '\t\tGap Extend = ', ge

  # Ensure scoring values are valid
  if ms < 0 or mms > 0 or go > 0 or ge > 0:
    print 'ERROR: Bad scoring parameter values'
    sys.exit(0)
  if ms + 3 * mms > 0:
    print 'ERROR: Positive expected score for random DNA sequences'
    print ms - 3*mms
    sys.exit(0)


  return locAL(seq1, seq2, ms, mms, go, ge, silenced = False, bestseq1 = True)

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

def external(seq1, seq2, mscore, mmscore, goscore, gescore):
  # Used when calling locAL from another python script
  return locAL(seq1, seq2, mscore, mmscore, goscore, gescore, silenced = True, bestseq1 = False)

def external_bestseq1(seq1, seq2, mscore, mmscore, goscore, gescore):
  # Finds the best alignment that includes all of seq1
  return locAL(seq1, seq2, mscore, mmscore, goscore, gescore, silenced = False, bestseq1 = True)

def locAL(seq1, seq2, ms, mms, go, ge, silenced = True, bestseq1 = False):
  # bestseq1 will force all of seq1 to be in the alignment
  s_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]
  d_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]
  bt_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]

  # initialize d and i, which can have negative scores
  d_mat[1][0] = go + ge
  d_mat[0][1] = go + ge
  for i in xrange(2, len(seq1) + 1):
    d_mat[i][0] = d_mat[i - 1][0] + ge
  for j in xrange(2, len(seq2) + 1):
    d_mat[0][j] = d_mat[0][j - 1] + ge
  i_mat = copy.deepcopy(d_mat)

  if bestseq1:
    s_mat[1][0] = go + ge
    for i in range(2, len(seq1) + 1):
      s_mat[i][0] = s_mat[i-1][0] + ge

  best = (0, 0, 0)

  if not silenced:
    print 'Generating alignment... Sequence len:', len(seq1)
  for i in xrange(1, len(seq1) + 1):
    for j in xrange(1, len(seq2) + 1):
      diag = s_mat[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], ms, mms, go, ge)

      d_up = d_mat[i - 1][j] + ge
      s_up = s_mat[i - 1][j] + go + ge
      d_mat[i][j] = max(d_up, s_up)

      i_left = i_mat[i][j - 1] + ge
      s_left = s_mat[i][j - 1] + go + ge
      i_mat[i][j] = max(i_left, s_left)

      if not bestseq1:
        s_mat[i][j] = max(0, d_mat[i][j], i_mat[i][j], diag)
      if bestseq1:
        s_mat[i][j] = max(d_mat[i][j], i_mat[i][j], diag)

      # update bt matrix
      if not bestseq1:
        if max(0, d_mat[i][j], i_mat[i][j], diag) == 0:
          bt_mat[i][j] = bt_zero
        if max(0, d_mat[i][j], i_mat[i][j], diag) == d_mat[i][j]:
          bt_mat[i][j] = bt_up
        if max(0, d_mat[i][j], i_mat[i][j], diag) == i_mat[i][j]:
          bt_mat[i][j] = bt_left
        if max(0, d_mat[i][j], i_mat[i][j], diag) == diag:
          bt_mat[i][j] = bt_diag
      if bestseq1:
        if max(d_mat[i][j], i_mat[i][j], diag) == d_mat[i][j]:
          bt_mat[i][j] = bt_up
        if max(d_mat[i][j], i_mat[i][j], diag) == i_mat[i][j]:
          bt_mat[i][j] = bt_left
        if max(d_mat[i][j], i_mat[i][j], diag) == diag:
          bt_mat[i][j] = bt_diag

      if not bestseq1:
        if s_mat[i][j] >= best[0]:
          best = (s_mat[i][j], i, j)

  if bestseq1:
    for j in range(1, len(seq2) + 1):
      if s_mat[len(seq1)][j] >= best[0]:
        best = (s_mat[i][j], len(seq1), j)

    for i in range(1, len(seq1) + 1):
      bt_mat[i][0] = bt_up

  (alignLen, matches, mismatches, numgaps, numGapExtends) = printAlignment(seq1, seq2, bt_mat, best[1:], silenced)
  if not silenced: 
    print 'Score =', best[0]
    if best[0] == 0:
      print 'SUGGESTION: Use different scoring parameters. No positive scoring alignment was found.'
  # print s_mat, '\n'
  # print bt_mat, '\n'
  # print d_mat, '\n'
  # print i_mat

  return (alignLen, matches, mismatches, numgaps, numGapExtends, best[1:])

def printAlignment(seq1, seq2, btMatrix, best_xy, silenced = True):
  # Given the backtracking matrix, prints out the best local alignment
  (i, j) = best_xy
  query = ''
  db = ''
  numgaps = 0
  mismatches = 0
  matches = 0

  while btMatrix[i][j] != bt_zero:
    if btMatrix[i][j] == bt_diag:
      query = string.join((query, seq1[i - 1]), '')
      db = string.join((db, seq2[j - 1]), '')
      if seq1[i - 1] == seq2[j - 1]:
        matches += 1
      else:
        mismatches += 1
      i, j = i - 1, j - 1

    if btMatrix[i][j] == bt_up:
      query = string.join((query, seq1[i - 1]), '')
      db = string.join((db, '-'), '')
      numgaps += 1
      i, j = i - 1, j

    if btMatrix[i][j] == bt_left:
      query = string.join((query, '-'), '')
      db = string.join((db, seq2[j - 1]), '')
      numgaps += 1
      i, j = i, j - 1      

  query = query[::-1]
  db = db[::-1]

  (dels, ins) = printSeqs(query, db, silenced)
  totalnumgaps = dels[0] + ins[0]
  totalgap_positions = dels[1] + ins[1]
  numGapExtends = totalgap_positions - totalnumgaps
  avg_gap_len = 0
  avg_ins_len = 0
  avg_del_len = 0
  if totalnumgaps > 0:
    avg_gap_len = float(totalgap_positions) / float(totalnumgaps)
  if ins[0] > 0:
    avg_ins_len = float(ins[1]) / float(ins[0])
  if dels[0] > 0:
    avg_del_len = float(dels[1])/ float(dels[0])

  if not silenced:
    print 'Alignment Len:', len(query), '\nMatches:', matches, '\tMismatches:', mismatches
    print 'Gaps:', totalnumgaps, '\tAvg Gap Len:', avg_gap_len
    print 'Ins:', ins[0], ' \tAvg Ins Len:', avg_ins_len
    print 'Del:', dels[0], ' \tAvg Del Len:', avg_del_len

  return (len(query), matches, mismatches, numgaps, numGapExtends)

def printSeqs(query, db, silenced = True):
  # Returns information about gaps
  m = ''
  dels = [0, 0]   # starts, total
  ins = [0, 0]    # starts, total
  numGapExtends = 0
  for i in range(len(query)):
    if query[i] == db[i]:
      m += '|'
    else:
      if query[i] == '-':
        m += ' '
        dels[1] += 1
        if query[i-1] != '-':
          dels[0] += 1
      if db[i] == '-':
        m += ' '
        ins[1] += 1
        if db[i-1] != '-':
          ins[0] += 1
    if query[i] != db[i] and query[i] != '-' and db[i] != '-':
      m += '*'
  if not silenced:
    print query, '\n', m, '\n', db
  return (dels, ins)

def score(base1, base2, ms, mms, go, ge):
  # Returns the numeric alignment score between the input bases
  if base1 == 'go' or base2 == 'go':
    return go
  if base1 == 'ge' or base2 == 'ge':
    return ge
  if base1 == base2:
    return ms
  else:
    return mms


if __name__ == '__main__':
  # main
  start = datetime.datetime.now()
  print 'Start:', start
  main()
  end = datetime.datetime.now()
  print '\nEnd:', end, '\nTotal:', end - start
