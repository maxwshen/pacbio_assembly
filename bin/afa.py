# Alignment-Free Assembly by finding the most read-supported path

import sys, string, datetime, random, copy, os, commands, fnmatch
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

  ktmer = 'TTTTCATCTGGAAGCGATTCTCGGGGA'   # pos = 4478401
  print gr.find(ktmer)

  traversed = TraversedReads(ktmer, creads, headers)
  while True:
    ktmer = move(ktmer, traversed, creads, headers, gr)
    print gr.find(ktmer)

  return

def move(ktmer, traversed, creads, headers, gr):
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

  def get_best_candidate(cands, traversed, creads, headers):
    scores = []
    for ckmer in cands:
      for r in headers[ckmer]:
        score = 0
        if r in traversed.tr:
          score += (get_pos_in_read(ckmer, creads[r]) - traversed.get(r)) ** 2
      scores.append(score)
    sorted_cands = [x for (y,x) in sorted(zip(scores, cands))]
    print '  ', sorted(scores), '\n  ', [gr.find(s) for s in sorted_cands]
    return sorted_cands[-1]

  cands, dists = find_candidates(ktmer, creads, headers)
  # curr_reads = set(headers[ktmer])
  # cands.sort(key = lambda x : shared_reads(curr_reads, headers[x]), reverse = True)
  if len(cands) > 0:
    new_ktmer = get_best_candidate(cands, traversed, creads, headers)
  else:
    new_ktmer = ''
    print 'No candidates found'

  traversed.update(new_ktmer, creads, headers)
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
