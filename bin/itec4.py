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
  ec_tool = '/home/mshen/research/bin/error_correction_1218.sh'
  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', ktmer_headers_file, '\nEC Tool:', ec_tool

  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool)


def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool):
  overlap_accuracy_cutoff = 80
  overlap_length_cutoff = 200
  num_attempts = 10
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  contigs = []
  curr_ktmer = ktmers[0]
  ktmers = ktmers[1:]
  h = get_read_with_most_neighbors(curr_ktmer, headers, creads)
  curr_contig = [error_correct(ec_tool, h, headers, creads, hr, rr)]
  traversed_headers = [h]

  # MAIN LOOP
  for direction in ['right', 'left']:
    counter = 0
    while True:
      counter += 1
      print 'iteration', counter, direction
      old_h = h
      stop = False
      for i in range(num_attempts):
        print 'Attempt', i
        possible_heads = extend_n(h, headers, creads, traversed_headers, direction)
        print len(possible_heads), 'candidates for extension'
        if len(possible_heads) == 0:
          print 'could not extend further'
          continue
        for head in possible_heads:
          candidate_read = rr[hr.index(head)]
          if test_overlap(candidate_read, curr_contig[-1], overlap_accuracy_cutoff, overlap_length_cutoff):
            h = head
            consensus_temp = error_correct(ec_tool, head, headers, creads, hr, rr)
            if len(consensus_temp) == 0:
              print 'failed to error correct'
              continue
            if direction == 'right':
              curr_contig.append(consensus_temp)
            if direction == 'left':
              curr_contig.insert(0, consensus_temp)
            traversed_headers.append(h)
            break
        if h == old_h:
          for p in possible_heads:
            traversed_headers.append(p)
        else:
          break
      if h == old_h:
        print 'No new reads overlapped with current contig'
        break

  # ASSESS RESULTS
  print old_h
  contig = ''
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'   # Concise output
  for k in curr_contig:
    if k != '':
      contig += '>c1\n' + k + '\n'
  with open('out.fasta', 'w') as f:
    f.write(contig)
  status = commands.getstatusoutput(blasr_exe + ' out.fasta ' + e_coli_genome + ' ' + blasr_options + ' > out.txt')[1]


def test_overlap(seq1, seq2, acc_cutoff, len_cutoff):
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'   # Concise output

  temps1 = 'temp_seq1.fasta'
  temps2 = 'temp_seq2.fasta'
  with open(temps1, 'w') as f:
    f.write('>1\n' + seq1)
  with open(temps2, 'w') as f:
    f.write('>2\n' + seq2)

  status = commands.getstatusoutput(blasr_exe + ' ' + temps1 + ' ' + temps2 + ' ' + blasr_options)[1]
  if len(status.strip()) == 0:
    return False
  print status
  accuracy = float(status.split()[5])
  length = int(status.split()[-1])
  end_align_r1 = int(status.split()[7])
  total_len_r1 = int(status.split()[8])
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_pos_r2 = int(status.split()[9])
  print 'accuracy:', accuracy
  dist_from_end = 200
  return accuracy >= acc_cutoff and length > len_cutoff and beg_pos_r2 < dist_from_end and end_pos_r1 < dist_from_end


def extend_n(header, headers, creads, traversed_headers, direction):
  dist_to_end = get_dist_to_end(header, creads, direction)
  ktmers = dist_to_end.keys()
  print len(ktmers), 'ktmers in 1-deg nhood'

  reads = []
  num_neighbors = []
  for k in ktmers:
    next_read = find_extending_read(k, headers, creads, dist_to_end[k], direction)
    if len(next_read) != 0:
      accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
      reads += accepted
      for nr in next_read:
        num_neighbors.append(len(creads[nr]) / 2 - 1)
  if len(reads) != 0:
    return [x for (y, x) in sorted(zip(num_neighbors, reads), reverse = True)]

  # Try finding a read that comes close to passing the current read
  for k in dist_to_end.keys():
    dist_to_end[k] -= 100

  print 'trying n-degree nhood extension'
  traversed_ktmers = set(ktmers)
  while True:
    print '\tkmers considered:', len(traversed_ktmers)
    new_ktmers = defaultdict(list)   # Key = new-ktmer, Val = [old ktmer that is connected]
    for kt in ktmers:
      for kn in find_neighboring_ktmers(kt, headers, creads):
        if kn not in traversed_ktmers:
          if kn not in new_ktmers.keys() or kt not in new_ktmers[kn]:
            new_ktmers[kn].append(kt)
    if len(new_ktmers.keys()) == 0:
      print 'None found'
      return ''

    print 'found', len(new_ktmers.keys()), 'ktmers from', len(ktmers), 'ktmers'
    new_dist_to_end = dict()
    reads = []
    num_neighbors = []
    for kt in new_ktmers.keys():
      k = new_ktmers[kt][0]
      dist = dist_bw_ktmers(kt, k, headers, creads)
      if dist == None:
        continue
      if direction == 'right':
        new_dist_to_end[kt] = dist_to_end[k] + dist
      if direction == 'left':
        new_dist_to_end[kt] = dist_to_end[k] - dist
      next_read = find_extending_read(kt, headers, creads, new_dist_to_end[kt], direction)
      if len(next_read) != 0:
        accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
        reads += accepted
        for nr in next_read:
          num_neighbors.append(len(creads[nr]) / 2 - 1)
    if len(reads) != 0:
      return [x for (y, x) in sorted(zip(num_neighbors, reads), reverse = True)]

    # Filter out those who are too far
    for k in ktmers:
      traversed_ktmers.add(k)

    # print new_dist_to_end
    ktmers = []
    limit = 1000       # Don't backtrack too far
    for kt in new_ktmers.keys():
      if new_dist_to_end[kt] < limit and kt not in traversed_ktmers:
        ktmers.append(kt)

    dist_to_end = new_dist_to_end


def find_neighboring_ktmers(ktmer, headers, creads):
  # Current speed bottleneck 12/25/14
  neighbors = []
  for h in headers[ktmer]:
    for i in range(len(creads[h])):
      if i % 2 == 1:
        neighbors.append(creads[h][i])
  return neighbors


def get_dist_to_end(header, creads, direction):
  dist_to_end = dict()    # Key = ktmer in read, Val = distance to end of read
  if direction == 'right':
    for i in range(len(creads[header])):
      if i % 2 == 1:
        dist_to_end[creads[header][i]] = 0
      elif i > 0:
        for k in dist_to_end:
          dist_to_end[k] += int(creads[header][i])
  if direction == 'left':
    curr = 0
    for i in range(len(creads[header])):
      if i % 2 == 1:
        dist_to_end[creads[header][i]] = curr
      else:
        curr += int(creads[header][i])
  return dist_to_end


def dist_bw_ktmers(kt1, kt2, headers, creads):
  # Finds distance between kt1 and kt2 in a read. If kt2 is before kt1, dist is negative
  found = False
  for h in headers[kt1]:
    if kt2 in creads[h]:
      header = h
      found = True

  if not found:
    return None

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


def find_extending_read(ktmer, headers, creads, dist, direction):
  # If a read in ktmer extends past dist, return it
  valid = []
  if direction == 'right':
    for h in headers[ktmer]:
      curr_dist = 0
      track = False
      for i in range(len(creads[h])):
        if creads[h][i] == ktmer:
          track = True
        if track and i % 2 == 0:
          curr_dist += int(creads[h][i])
      if curr_dist > dist:
        valid.append(h)
  if direction == 'left':
    for h in headers[ktmer]:
      curr_dist = 0
      track = True
      for i in range(len(creads[h])):
        if creads[h][i] == ktmer:
          track = False
        if track and i % 2 == 0:
          curr_dist += int(creads[h][i])
      if curr_dist > dist:
        valid.append(h)
  return valid


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

  reads = [header, rr[hr.index(header)]]
  for ch in collected_h:
    reads.append(ch)
    reads.append(rr[hr.index(ch)])
  print len(reads) / 2, 'reads used for correction'

  temp_orig_file = 'temp_orig.fasta'
  with open(temp_orig_file, 'w') as f:
    f.write(header + '\n' + rr[hr.index(header)])

  temp_nhood_file = 'temp_nhood.fasta'
  with open(temp_nhood_file, 'w') as f:
    f.write('\n'.join(reads))

  ec_out = 'corrected_' + temp_orig_file
  status = commands.getstatusoutput(ec_tool + ' ' + temp_orig_file + ' ' + temp_nhood_file)[1]
  if 'ERROR' in status:
    print status
    return ''

  with open(ec_out, 'r') as f:  
    consensus = f.readlines()[1].strip()
  print 'consensus len:', len(consensus), 'out of', len(rr[hr.index(header)])
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