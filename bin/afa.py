# Alignment-Free Assembly by finding the most read-supported path

import sys, string, datetime, random, copy, os, commands, fnmatch, re
import numpy as np
from collections import defaultdict

import mylib as ml
import find_read
import convert_creads_to_nhoods
import itec4

prior = '/home/mshen/research/'
reads_fn = prior + 'data/reads.20k.rc.fasta'
creads_fn = prior + 'data/20k_v2/temp_creads.outrx_27_6_rc_v2.out'
ktmer_headers_fn = prior + 'data/20k_v2/temp_ktmer_headersrx_27_6_rc_v2.out'
genome_fn = prior + 'data/ecoli_consensus_mark.fasta'

bt_fraction = 0.01    # If the current best score is less than fraction * best last score, backtrack
forget_cutoff = 50000
traversed = []
score_history = ([0], [0])    # First is score history, second is num_candidates history

hr, rr = ml.read_fasta(reads_fn)
headers = itec4.build_headers_dict(ktmer_headers_fn)
creads = itec4.build_creads_dict(creads_fn, hr, rr)
for i in range(len(hr)):
  hr[i] = hr[i].split()[0]

def main():
  print 'Reads File:', reads_fn, '\ncreads File:', creads_fn, '\nktmer Headers File:', \
    ktmer_headers_fn

  afa(reads_fn, ktmer_headers_fn, creads_fn)
  return

def afa(reads_fn, ktmer_headers_fn, creads_fn):
  gh, gr = ml.read_fasta(genome_fn)
  gr = gr[0]

  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  num_shared_cutoff = 4
  # ktmer = 'TTTTCATCTGGAAGCGATTCTCGGGGA'   # pos = 4478401
  # ktmer = 'TGGTGTTTTTGAGGTTCTCCAGTGGCT'   # pos = 4506698
  # ktmer = 'CGCCGCCAATCGCCAGACCCAGACGGC'     # pos = 4581355 -> 100k, skips over backtracking part
  ktmer = 'CAAAGGAGTCCAGGCTTCCCCCATAAC'     # pos = 68407, repeat issue at 101k

  # ktmer = 'AACCGCCTTGTCCGAGAGTAAAGCTTG'
  # ktmer = 'TCGTACATTTAAGAAATTAAATCATTT'
  # ktmer = 'TGTTCTGTAACGCCAGCGAATCGGTAT'

  find_ktmer_in_genome(ktmer, gr)

  seenReads = TraversedReads()
  seenReads.update(ktmer)
  traversed.append(ktmer)
  check_bt = False
  while True:
    ktmer = move(ktmer, seenReads, gr, num_shared_cutoff)
    if ktmer == 'BACKTRACK':
      print 'backtracking 3...'
      ktmer = traversed[-2]
      if check_bt:
        print 'Backtracking did not help'
        break
      check_bt = True
    else:
      check_bt = False
    find_ktmer_in_genome(ktmer, gr)
    print '  curr pos in reads:', [get_pos_in_read(ktmer, creads[r]) for r in headers[ktmer]]
    print '  len reads:', [get_len_cread(creads[r]) for r in headers[ktmer]]

  return

def check_pct_genomic(ktmers, genome):
  found = 0
  for i in range(len(ktmers)):
    kt = ktmers[i]
    if genome.find(kt) != -1 and genome.find(ml.reverse_complement(kt)):
      found += 1
    print i, found, float(found) / float(i + 1)
  print '# Genomic:', found
  print 'pct Genomic:', float(found) / float(len(ktmers))

def find_ktmer_in_genome(ktmer, genome):
  print ', '.join([str(s.start()) for s in re.finditer(ktmer, genome)]), ktmer
  print ', '.join([str(s.start()) for s in re.finditer(ml.reverse_complement(ktmer), genome)]), ml.reverse_complement(ktmer)
  return

def find_repeats(ktmer):
  def find_index_after_lastrepeat(cread):
    ktm = dict()
    for i in range(len(cread)):
      if i % 2 == 1:
        kt = cread[i]
        if kt not in ktm:
          ktm[kt] = 1
        else:
          ktm[kt] += 1
    last_pos = -1
    for i in range(len(cread) - 1, 0, -1):
      if i % 2 == 1:
        kt = cread[i]
        if ktm[kt] > 1:
          last_pos = i + 2
          break
    return last_pos

  # BEGIN def find_repeats(...)
  next_ktmers = []
  next_dists = []
  for r in headers[ktmer]:
    cc = creads[r]
    instances = [i for i, x in enumerate(cc) if x == ktmer]
    if len(instances) > 1:
      nkti = find_index_after_lastrepeat(cc)
      if nkti < len(cc) and cc[nkti] not in next_ktmers:
        next_ktmers.append(cc[nkti])
        next_dists.append(get_dist_in_cread(cc, cc.index(ktmer), nkti))
  return next_ktmers, next_dists

def move(ktmer, seenReads, gr, num_shared_cutoff):
  def find_candidates(ktmer):
    next_ktmers, next_dists = find_repeats(ktmer)
    if len(next_ktmers) > 0:
      print 'Found repeat kt-mer'
      return next_ktmers, next_dists
    else:
      # If this kt-mer is not repeated in reads
      reads = headers[ktmer]
      next_ktmers = []
      next_dists = []
      for r in reads:
        nkti = creads[r].index(ktmer) + 2
        if nkti < len(creads[r]) and creads[r][nkti] not in next_ktmers:
          next_ktmers.append(creads[r][nkti])
          next_dists.append(int(creads[r][nkti - 1]))
      return next_ktmers, next_dists

  def shared_reads(base_reads, new_reads):
    return len(base_reads.intersection(new_reads))

  def score(kmer_pos, first_pos):
    if kmer_pos > first_pos:
      return (kmer_pos - first_pos) ** 2
    else:
      return 0

  def get_best_candidate(cands, seenReads, num_shared_cutoff):
    if len(cands) == 0:
      print 'No candidates found'
      return ''

    # # If one-overlap, immediately take it
    # for i in range(len(cands)):
    #   if cands[i][:-1] == ktmer[1:]:
    #     return cands[i]

    scores = []
    for ckmer in cands:
      scoresum = 0
      num_shared = 0
      for r in headers[ckmer]:
        if r in seenReads.tr:
          num_shared += 1
          scoresum += score(get_pos_in_read(ckmer, creads[r]), seenReads.get(r))
      if num_shared <= num_shared_cutoff:
        scoresum = 0
      scores.append(scoresum)
    if sum(scores) == 0:
      return ''

    sorted_cands = [x for (y,x) in sorted(zip(scores, cands))]
    print '  cands scores:', sorted(scores), '\n  cands in genome:', [gr.find(s) for s in sorted_cands]
    # print '  ', sorted_cands

    best_score = sorted(scores)[-1]
    if score_history[1][-1] > 1 and best_score < score_history[0][-1] * bt_fraction:
      print 'backtracking...'
      return 'BACKTRACK'
    else:
      score_history[0].append(best_score)
      score_history[1].append(len(sorted_cands))
      return sorted_cands[-1]

  # BEGIN def move(...)
  cands, dists = find_candidates(ktmer)
  new_cands = [s for s in cands if s not in traversed]
  new_ktmer = get_best_candidate(new_cands, seenReads, num_shared_cutoff)
  if new_ktmer == 'BACKTRACK':
    print 'backtracking 2...'
    return 'BACKTRACK'
  while new_ktmer == '':
    num_shared_cutoff -= 1
    if num_shared_cutoff < 0:
      print 'num_shared_cutoff is less than 0 and still no candidates were found'
      sys.exit(0)
    new_ktmer = get_best_candidate(new_cands, seenReads, num_shared_cutoff)
    if new_ktmer == 'BACKTRACK':
      print 'backtracking 2.5...'
      return 'BACKTRACK'

  new_dist = dists[cands.index(new_ktmer)]
  seenReads.update(new_ktmer)
  nf = seenReads.forget(new_dist, forget_cutoff)
  print 'Moved forward by', new_dist, 'and forgot', nf, 'reads'
  traversed.append(new_ktmer)
  return new_ktmer

def get_pos_in_read(kmer, cread):
  if kmer in cread:
    return sum([int(cread[s]) for s in range(cread.index(kmer)) if s % 2 == 0])
  else:
    print 'error:', kmer, 'not in', cread

def get_len_cread(cread):
  return sum([int(cread[s]) for s in range(len(cread)) if s % 2 == 0])

def get_dist_in_cread(cread, index1, index2):
  return sum([int(cread[s]) for s in range(index1, index2 + 1) if s % 2 == 0  ])

class TraversedReads():
  def __init__(self):
    self.tr = dict()
    self.td = dict()

  def forget(self, dist, cutoff):
    num_forgotten = 0
    for k in self.td.keys():
      self.td[k] += dist
      if self.td[k] > cutoff:
        del self.td[k]
        del self.tr[k]
        num_forgotten += 1
    return num_forgotten

  def update(self, kmer):
    for read in headers[kmer]:
      pos = get_pos_in_read(kmer, creads[read])
      if read not in self.tr:
        self.tr[read] = pos
        self.td[read] = 0
      else:
        if pos < self.tr[read]:
          self.tr[read] = pos
    return

  def get(self, read):
    if read in self.tr:
      return self.tr[read]
    else:
      print read, 'not in TraversedReads'
      return -1

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
