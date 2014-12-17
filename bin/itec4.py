# Treats reads as the basis for error correction and extension, rather than kt-mers as in iterative_ec3.py
#
# Related files:
# combine_ec_contigs.py : Combines the read sets produced by this code into contigs

import sys, string, datetime, random, copy, os, commands
import numpy as np
from collections import defaultdict

import read_fasta as rf
import find_read

def main():
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  creads_file = '/home/mshen/research/data/22.4_creads.out'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  ec_tool = '/home/mshen/research/bin/consensus_correction.sh'
  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', ktmer_headers_file, '\nEC Tool:', ec_tool

  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool)

def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool):
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  contigs = []
  curr_contig = []

  curr_ktmer = ktmers[0]
  ktmers = ktmers[1:]
  h = get_read_with_most_neighbors(curr_ktmer, headers, creads)
  traversed_headers = [h]
  # error_correct(ec_tool, h, headers, creads, hr, rr)
  i = 0
  while len(h) != 0:
    i += 1
    print i
    old_h = h
    h = extend_right(h, headers, creads, traversed_headers)
    if len(h) == 0:
      h = extend_right_2(old_h, headers, creads, traversed_headers)
    traversed_headers.append(h)

  contig = ''
  for k in traversed_headers:
    if k != '':
      contig += k + '\n' + rr[hr.index(k)] + '\n'
  with open('out.fasta', 'w') as f:
    f.write(contig)

def extend_right(header, headers, creads, traversed_headers):
  dist_to_end = dict()    # Key = ktmer in read, Val = distance to end of read
  for i in range(len(creads[header])):
    if i % 2 == 1:
      dist_to_end[creads[header][i]] = 0
    elif i > 0:
      for k in dist_to_end:
        dist_to_end[k] += int(creads[header][i])
  ktmers = dist_to_end.keys()
  print len(ktmers), 'ktmers'
  print dist_to_end

  for k in ktmers:
    next_read = find_extending_read_rightend(k, headers, creads, dist_to_end[k])
    if len(next_read) != 0 and next_read not in traversed_headers:
      return next_read

  print 'not found'
  return ''

def extend_right_2(header, headers, creads, traversed_headers):
  dist_to_end = dict()    # Key = ktmer in read, Val = distance to end of read
  for i in range(len(creads[header])):
    if i % 2 == 1:
      dist_to_end[creads[header][i]] = 0
    elif i > 0:
      for k in dist_to_end:
        dist_to_end[k] += int(creads[header][i])
  ktmers = dist_to_end.keys()
  print len(ktmers), 'ktmers'
  print dist_to_end

  for k in ktmers:
    for h in headers[k]:
      dist_to_k = dict()    # Key = ktmer in read, Val = distance to k
      kts = []
      for i in range(len(creads[h])):
        if i % 2 == 1:
          kts.append(creads[h][i])
      kts.remove(k)
      for kt in kts:
        dist_to_k[kt] = dist_bw_ktmers(kt, k, h, creads)
      for kt in dist_to_k.keys():
        next_read = find_extending_read_rightend(kt, headers, creads, dist_to_end[k] + dist_to_k[kt])
        if len(next_read) != 0 and next_read not in traversed_headers:
          return next_read

  print 'not found'
  return ''

def dist_bw_ktmers(kt1, kt2, header, creads):
  # Finds distance between kt1 and kt2 in a read. If kt2 is before kt1, dist is negative
  dist = 0
  grab = False
  for i in range(len(creads[header])):
    if creads[header][i] == kt1 or creads[header][i] == kt2:
      if grab == True:
        grab = False
        break
      if grab == False:
        grab = True
    if grab and i % 2 == 0:
      dist += int(creads[header][i])
  if creads[header].index(kt1) > creads[header].index(kt2):
    dist *= -1
  return dist

def find_extending_read_rightend(ktmer, headers, creads, dist):
  # If a read in ktmer extends past dist, return it
  for h in headers[ktmer]:
    curr_dist = 0
    track = False
    for i in range(len(creads[h])):
      if creads[h][i] == ktmer:
        track = True
      if track and i % 2 == 0:
        curr_dist += int(creads[h][i])
    if curr_dist > dist:
      return h
  return ''

def error_correct(ec_tool, header, headers, creads, hr, rr):
  reads = []
  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    return ''
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      collected_h.add(h)
  print len(collected_h)

  reads = [header, rr[hr.index(header)]]
  for ch in collected_h:
    reads.append(ch)
    reads.append(rr[hr.index(ch)])
  print len(reads)

  with open('temp_orig.fasta', 'w') as f:
    f.write(header + '\n' + rr[hr.index(header)])

  temp_nhood_file = 'temp_nhood.fasta'
  with open(temp_nhood_file, 'w') as f:
    f.write('\n'.join(reads))

  ec_out = temp_nhood_file + '_consensus.fasta'
  status = commands.getstatusoutput(ec_tool + ' ' + temp_nhood_file)[1]

  with open(ec_out) as f:  
    consensus = f.readlines()[1].strip()
  print len(rr[hr.index(header)]), len(consensus)
  return consensus

def get_read_with_most_neighbors(ktmer, headers, creads):
  best_h = ''
  best_neighbors = 0
  for h in headers[ktmer]:
    num_neighbors = len(creads[h]) / 2 - 1
    print h, num_neighbors
    if num_neighbors > best_neighbors:
      best_neighbors = num_neighbors
      best_h = h
  return best_h

def build_creads_dict(creads_file, reads_file):
  creads = defaultdict(list)   # Key = header, Val = creads 
  h, r = rf.read_fasta(reads_file)
  with open(creads_file) as f:
    for i, line in enumerate(f):
      creads[h[i]] = line.split()
  return creads

def build_headers_dict(ktmer_headers_file):
  headers = defaultdict(list)   # Key = ktmer, Val = [headers]
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      headers[words[0]] = words[1:]
  return headers

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start