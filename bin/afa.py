# Alignment-Free Assembly by finding the most read-supported path

import sys, string, datetime, random, copy, os, commands, fnmatch, re
import numpy as np
from collections import defaultdict

import read_fasta as rf
import find_read
import convert_creads_to_nhoods
import itec4

prior = '/home/mshen/research/'
genome_fn = prior + 'data/ecoli_consensus_mark.fasta'

def main():
 
  reads_fn = prior + 'data/reads.20k.rc.fasta'
  creads_fn = prior + 'data/20k_v2/temp_creads.outrx_27_6_rc_v2.out'
  ktmer_headers_fn = prior + 'data/20k_v2/temp_ktmer_headersrx_27_6_rc_v2.out'

  print 'Reads File:', reads_fn, '\ncreads File:', creads_fn, '\nktmer Headers File:', \
    ktmer_headers_fn

  afa(reads_fn, ktmer_headers_fn, creads_fn)
  return

def afa(reads_fn, ktmer_headers_fn, creads_fn):
  hr, rr = rf.read_fasta(reads_fn)
  headers = itec4.build_headers_dict(ktmer_headers_fn)
  creads = itec4.build_creads_dict(creads_fn, hr, rr)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  gh, gr = rf.read_fasta(genome_fn)
  gr = gr[0]

  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  num_shared_cutoff = 6
  ktmer = 'TTTTCATCTGGAAGCGATTCTCGGGGA'   # pos = 4478401
  print gr.find(ktmer)

  seenReads = TraversedReads(ktmer, creads, headers)
  traversed = set(ktmer)
  while True:
    ktmer = move(ktmer, seenReads, traversed, creads, headers, gr, num_shared_cutoff)
    # print gr.find(ktmer)
    print 'pos: ', ', '.join([str(s.start()) for s in re.finditer(ktmer, gr)])

  return

def move(ktmer, seenReads, traversed, creads, headers, gr, num_shared_cutoff):
  def find_candidates(ktmer, creads, headers):
    reads = headers[ktmer]
    next_ktmers = []
    next_dists = []
    for r in reads:
      nkti = creads[r].index(ktmer) + 2
      if nkti < len(creads[r]) and creads[r][nkti] not in next_ktmers:
        next_ktmers.append(creads[r][nkti])
        next_dists.append(creads[r][nkti - 1])
    return next_ktmers, next_dists

  def shared_reads(base_reads, new_reads):
    return len(base_reads.intersection(new_reads))

  def score(kmer_pos, first_pos):
    if kmer_pos > first_pos:
      return (kmer_pos - first_pos) ** 2
    else:
      return 0

  def get_best_candidate(cands, seenReads, creads, headers, num_shared_cutoff):
    if len(cands) == 0:
      print 'No candidates found'
      return ''

    scores = []
    for ckmer in cands:
      scoresum = 0
      num_shared = 0
      for r in headers[ckmer]:
        if r in seenReads.tr:
          num_shared += 1
          scoresum += score(get_pos_in_read(ckmer, creads[r], seenReads.get(r)))
      if num_shared <= num_shared_cutoff:
        scoresum = 0
      scores.append(scoresum)
    sorted_cands = [x for (y,x) in sorted(zip(scores, cands))]
    print '  ', sorted(scores), '\n  ', [gr.find(s) for s in sorted_cands]
    return sorted_cands[-1]


  cands, dists = find_candidates(ktmer, creads, headers)
  cands = [s for s in cands if s not in traversed]
  new_ktmer = get_best_candidate(cands, seenReads, creads, headers, num_shared_cutoff)
  while new_ktmer == '':
    num_shared_cutoff -= 1
    if num_shared_cutoff < 0:
      print 'num_shared_cutoff is less than 0 and still no candidates were found'
      sys.exit(0)
    new_ktmer = get_best_candidate(cands, seenReads, creads, headers, num_shared_cutoff)

  seenReads.update(new_ktmer, creads, headers)
  traversed.add(new_ktmer)
  return new_ktmer

def get_pos_in_read(kmer, cread):
  if kmer in cread:
    return sum([int(cread[s]) for s in range(cread.index(kmer)) if s % 2 == 0])
  else:
    print 'error:', kmer, 'not in', cread


class TraversedReads():
  def __init__(self, kmer, creads, headers):
    self.tr = dict()
    self.update(kmer, headers, creads)

  def update(self, kmer, creads, headers):
    for read in headers[kmer]:
      pos = get_pos_in_read(kmer, creads[read])
      if read not in self.tr:
        self.tr[read] = pos
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
