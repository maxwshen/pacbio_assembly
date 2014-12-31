# Treats reads as the basis for error correction and extension, rather than kt-mers as in iterative_ec3.py
#
# Related files:
# combine_ec_contigs.py : Combines the read sets produced by this code into contigs

import sys, string, datetime, random, copy, os, commands
import numpy as np
from collections import defaultdict

import read_fasta as rf
import find_read
import kmer_matching

def main():
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  creads_file = '/home/mshen/research/data/22.4_creads.out'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  ec_tool = '/home/mshen/research/bin/error_correction_1218.sh'
  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', ktmer_headers_file, '\nEC Tool:', ec_tool

  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool)
  # ktmer_reads_pct_overlap(ktmer_headers_file, reads_file)


def ktmer_reads_pct_overlap(ktmer_headers_file, reads_file):
  # Finds read clusters (aligned to the genome) for all kt-mers
  def within(beg1, end1, beg2, end2):
    # Test if 2 is in 1
    if beg1 < beg2 < end1:
      return True
    if beg1 < end2 < end1:
      return True
    return False

  def expand(beg1, end1, beg2, end2):
    newbeg1 = beg1
    newend1 = end1
    if beg2 < beg1:
      newbeg1 = beg2
    if end2 > end1:
      newend1 = end2
    return [newbeg1, newend1]

  def within_all(clusters):
    for i in range(len(clusters)):
      for j in range(len(clusters)):
        if i != j:
          if within(clusters[i][0], clusters[i][1], clusters[j][0], clusters[j][1]):
            return i, j
    return None

  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'   # Concise output

  for kt in headers.keys():
    print kt
    clusters = []
    for h in headers[kt]:
      tempfile = 'temp.fasta'
      with open(tempfile, 'w') as f:
        f.write('>1\n' + rr[hr.index(h)]) 
      status = commands.getstatusoutput(blasr_exe + ' ' + tempfile + ' ' + e_coli_genome + ' ' + blasr_options)[1]
      if len(status) > 0:
        beg = int(status.split()[6])
        end = int(status.split()[7])
        found = False
        for c in clusters:
          if within(c[0], c[1], beg, end):
            found = True
            new = expand(c[0], c[1], beg, end)
            c[0] = new[0]
            c[1] = new[1]
            c[2] += 1
        if not found:
          clusters.append([beg, end, 1])
        # print clusters, beg, end

    # Do final clustering of clusters
    while True:
      found = False
      if within_all(clusters) is not None:
        found = True
        i, j = within_all(clusters)
        new = expand(clusters[i][0], clusters[i][1], clusters[j][0], clusters[j][1])
        clusters[i][0] = new[0]
        clusters[i][1] = new[1]
        clusters[i][2] += clusters[j][2]
        clusters.remove(clusters[j])
      if not found:
        break

    # Currently just prints genomic alignment of clusters.
    for c in clusters:
      print ' '.join([str(s) for s in c])

    # May want to write a sister function that doesn't rely on the genome


def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool):
  overlap_accuracy_cutoff = 75    # .
  overlap_length_cutoff = 200     # .
  num_attempts = 3                # Number of times to try nhood extension.
  limit_km_times_total = 1        # How many times to attempt k-mer matching extension per direction
  km_k = 15                       # .
  km_cutoff = 10                  # .
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  contigs = []

  num_contig_attempts = 400
  for i in range(num_contig_attempts):
    curr_ktmer = ktmers[i]

    h = get_read_with_most_neighbors(curr_ktmer, headers, creads)
    curr_contig = [error_correct(ec_tool, h, headers, creads, hr, rr)]
    curr_contig_headers = [h]
    traversed_headers = [h]

    # MAIN LOOP
    for direction in ['right', 'left']:
      counter = 0
      limit_km_times = limit_km_times_total
      while True:
        # Break condition: Current header doesn't change, meaning we couldn't find any extension candidates
        counter += 1
        print 'iteration', counter, direction
        old_h = h
        stop = False
        for i in range(num_attempts + limit_km_times):
          print 'Attempt', i
          km = False

          # Grab candidates via nhood extension or kmer matching
          if i < num_attempts:
            # Try nhood extension num_attempt times, then try kmer matching if still no extension
            possible_heads = extend_n(h, headers, creads, traversed_headers, direction)
            print len(possible_heads), 'candidates for extension'
            if len(possible_heads) == 0:
              print 'could not extend further'
              continue
          else:
            print 'Trying k-mer matching'
            limit_km_times -= 1
            km = True
            if direction == 'right':
              possible_heads = kmer_matching.kmer_matching(curr_contig[-1], reads_file, km_k, km_cutoff, file_bool = False)
            if direction == 'left':
              possible_heads = kmer_matching.kmer_matching(curr_contig[0], reads_file, km_k, km_cutoff, file_bool = False)
            if h in possible_heads:
              possible_heads.remove(h)

          # Test candidates for overlap
          for head in possible_heads:
            candidate_read = rr[hr.index(head)]
            if test_overlap(candidate_read, curr_contig[-1], overlap_accuracy_cutoff, overlap_length_cutoff, direction, relaxed = km):
              h = head
              consensus_temp = error_correct(ec_tool, head, headers, creads, hr, rr)
              if len(consensus_temp) == 0:
                print 'failed to error correct'
                continue
              if direction == 'right':
                curr_contig.append(consensus_temp)
                curr_contig_headers.append(h)
              if direction == 'left':
                curr_contig.insert(0, consensus_temp)
                curr_contig_headers.insert(0, h)
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
    contigs_fold = '/home/mshen/research/contigs/'
    contig_file = contigs_fold + 'contig_' + str(i) + '.fasta'
    contig_result = contigs_fold + 'contig_' + str(i) + 'results.fasta'
    for i in range(len(curr_contig)):
      if curr_contig[i] != '':
        contig += '>' + curr_contig_headers[i] + '\n' + curr_contig[i] + '\n'
    with open(contig_file, 'w') as f:
      f.write(contig)
    status = commands.getstatusoutput(blasr_exe + ' ' + contig_file +' ' + e_coli_genome + ' ' + blasr_options + ' > ' + contig_file)[1]


def test_overlap(base, candidate, acc_cutoff, len_cutoff, direction, relaxed = False):
  # Tests that seq1 is before seq2
  if direction == 'right':
    seq1 = base
    seq2 = candidate
  if direction == 'left':
    seq1 = candidate
    seq2 = base
  dist_from_end = 200
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
  beg_align_r1 = int(status.split()[6])
  end_align_r1 = int(status.split()[7])
  total_len_r1 = int(status.split()[8])
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[9])
  end_align_r2 = int(status.split()[10])
  total_len_r2 = int(status.split()[11])
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2          # Average alignment length
  print 'accuracy:', accuracy

  if not relaxed:
    if direction == 'right':
      return accuracy >= acc_cutoff and length > len_cutoff and end_pos_r1 < dist_from_end and end_pos_r2 > end_pos_r1
    if direction == 'left':
      return accuracy >= acc_cutoff and length > len_cutoff and beg_align_r2 < dist_from_end and beg_align_r1 > beg_align_r2
  else:
    return accuracy >= acc_cutoff and length > len_cutoff


def extend_n(header, headers, creads, traversed_headers, direction):
  leniency = 100    # If 1-degree nhood fails, we accept a read that extends beyond border - leniency
  backtrack_limit = 1000       # Don't backtrack too far
  num_kmers_cutoff = 500       # If we start considering this many kt-mers in nhood extension, stop.
  
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
    dist_to_end[k] -= leniency

  print 'trying n-degree nhood extension'
  traversed_ktmers = set(ktmers)
  while True:
    print '\tkmers considered:', len(traversed_ktmers)
    new_ktmers = defaultdict(list)   # Key = new-ktmer, Val = [old ktmer that is connected]
    if len(ktmers) > num_kmers_cutoff:
      print 'Stopping - too many kt-mers'
      return ''
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
    for kt in new_ktmers.keys():
      if new_dist_to_end[kt] < backtrack_limit and kt not in traversed_ktmers:
        ktmers.append(kt)

    dist_to_end = new_dist_to_end


def find_neighboring_ktmers(ktmer, headers, creads):
  # Speed bottleneck 12/25/14
  # But storing dict of all neighbors in memory is slower
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