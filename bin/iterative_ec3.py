# Performs iterative error correction and production of consensus sequences in contigs.
# Beginning with an arbitrary ktmer and new contig list,
#   1. Run error correction on this ktmer
#   2. Add consensus to the contig list
#   3. Find a new ktmer some distance away from the old ktmer and repeat
#     - If there are no new ktmers, end the contig and start a new contig list at a new ktmer

# Only uses the reads in (22,4)-mers as opposed to full neighborhoods, as in iterative_ec2.py
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
  ktmer_edges_file = '/home/mshen/research/data/22.4_ktmer_edges.out'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'

  reads_file = 'extracts_100k/extracted_reads_c2500000_s100000.fasta'
  ktmer_edges_file = '/home/mshen/research/temp_ktmer_edges2.out'
  ktmer_headers_file = '/home/mshen/research/temp_ktmer_headers2.out'
  min_dist = 20

  print 'Reads File:', reads_file, '\nktmer Edges File:', ktmer_edges_file, '\nktmer Headers File:', ktmer_headers_file 

  iterative_ec(reads_file, ktmer_headers_file, ktmer_edges_file, min_dist)

def put_ktmer_reads_into_file(ktmer, ktmer_headers_file, reads_file, out_file):
  headers = []
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      if line.split()[0] == ktmer:
        headers = line.split()[1:]
  with open(out_file, 'w+') as f:
    for i in range(len(headers)):
      if i % 3 == 0:
        f.write(find_read.find_read(headers[i], reads_file))
  return

def build_headers_dict(ktmer_headers_file):
  headers = defaultdict(list)   # Key = header, Val = [headers, index, dist_from_end, ...]
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      headers[words[0]] = words[1:]
  return headers

def get_all_reads_from_read(ktmer, edges, headers, reads_file, out_file):
  neighbors = [ktmer]
  for i in range(len(edges[ktmer])):
    if i % 2 == 0:
      neighbors.append(edges[ktmer][i])


  found_headers = []
  for kt in neighbors:
    heads = headers[kt]
    heads = [heads[i] for i in range(len(heads)) if i % 3 == 0]
    for h in heads:
      if h not in found_headers:
        found_headers.append(h)
  print found_headers

  with open(out_file, 'w') as f:
    f.write(find_read.find_reads(found_headers, reads_file))
  return

def iterative_ec(reads_file, ktmer_headers_file, ktmer_edges_file, min_dist):
  ec_tool = '/home/mshen/research/bin/consensus_correction.sh'

  ktmers = []
  with open(ktmer_edges_file) as f:
    for i, line in enumerate(f):
      ktmers.append(line.split()[0])
  random.shuffle(ktmers)

  print 'Building dicts...', datetime.datetime.now()
  edges = build_edges_dict(ktmer_edges_file)
  headers = build_headers_dict(ktmer_headers_file)
  print '...Done.', datetime.datetime.now()

  # temp_outfile = 'temp_allreadsfromread.fasta'
  # get_all_reads_from_read(ktmers[0], edges, headers, reads_file, temp_outfile)

  restart_contig = False
  contigs = []    # List of lists
  curr_contig = []
  curr = ktmers[0]
  direction = 'forward'
  orig = curr
  # while len(ktmers) > 0:
  while len(contigs) < 11:
    if restart_contig:
      print len(contigs), len(ktmers)
      restart_contig = False
      contigs.append(curr_contig)
      curr_contig = []
      curr = ktmers[0]
      direction = 'forward'
      orig = curr
    if curr in ktmers:
      del ktmers[ktmers.index(curr)]

    ktmer_reads_out = 'temp_ktmer_reads.fasta'
    put_ktmer_reads_into_file(curr, ktmer_headers_file, reads_file, ktmer_reads_out)

    ec_file = 'corrected_' + ktmer_reads_out
    status = commands.getstatusoutput(ec_tool + ' ' + ktmer_reads_out)[1]
    if 'ERROR' in status:
      print status
      restart_contig = True
      continue

    with open(ec_file) as f:
      ec_seq = f.readlines()[1].strip()
    if direction == 'forward':
      curr_contig.append(ec_seq)
    elif ec_seq not in curr_contig:   # Reinserts twice when we switch from forward to back
      curr_contig.insert(0, ec_seq)

    # max_dist = calc_max_dist(curr, ec_seq, direction)
    if direction == 'forward':
      max_dist = float('inf')
    else:
      max_dist = - float('inf')

    print curr, 'dist', max_dist, len(ec_seq), direction
    print edges[curr]

    found = False
    edge = find_best_edge(curr, edges, min_dist, max_dist, ktmers, headers, direction)
    # edge = find_random_edge(curr, edges, min_dist, max_dist, ktmers)
    if edge != '':
      curr = edge
    else:
      if direction == 'forward':
        direction = 'backward'
        curr = orig
      else:
        restart_contig = True

  contigs_file = 'ec_contigs.fasta'
  with open(contigs_file, 'w+') as f:
    for i in range(len(contigs)):
      for seq in contigs[i]:
        f.write('>' + str(i) + '\n' + seq + '\n')

  return

def calc_max_dist(curr, ec_seq, direction):
  if curr in ec_seq:
    max_dist = len(ec_seq) - ec_seq.index(curr)
    if direction == 'backward':
      max_dist -= len(ec_seq)
  else:
    max_dist = len(ec_seq) / 2
    if direction == 'backward':
      max_dist *= -1  
  return max_dist

def build_edges_dict(ktmer_edges_file):
  edges = defaultdict(list)
  with open(ktmer_edges_file) as f:
    for i, line in enumerate(f):
      edges[line.split()[0]] = line.split()[1:]
  return edges

def find_best_edge(ktmer, edges, mindist, max_dist, ktmers, headers, direction):
  # Returns a list of ktmers that are within min dist, max dist
  # direction = 'forward', or 'backward'
  if max_dist < 0:
    min_dist = mindist * -1
  else:
    min_dist = mindist
  curr_edges = edges[ktmer]
  best_edge = ''
  best_degree = 0
  num_candidates = 0
  num_passed_dirtest = 0
  for i in range(len(curr_edges)):
    if i % 2 == 1:
      if min_dist <= int(curr_edges[i]) <= max_dist or max_dist <= int(curr_edges[i]) <= min_dist:
        num_candidates += 1
        if direction_test(curr_edges[i - 1], headers, direction):
          num_passed_dirtest += 1
          if direction == 'forward':
            if len([s for s in edges[curr_edges[i - 1]] if s > min_dist]) > best_degree and curr_edges[i - 1] in ktmers:
              best_degree = len(edges[curr_edges[i - 1]])
              best_edge = curr_edges[i - 1]
          if direction == 'backward':
            if len([s for s in edges[curr_edges[i - 1]] if s < min_dist]) > best_degree and curr_edges[i - 1] in ktmers:
              best_degree = len(edges[curr_edges[i - 1]])
              best_edge = curr_edges[i - 1]


  print best_edge, best_degree / 2, headers[best_edge], 'candidates:', num_candidates, num_passed_dirtest
  return best_edge

def direction_test(ktmer, headers, direction):
  # Currently doesn't filter out reads on the wrong side. Naive - test if good enough
  dist_threshold = 500
  min_pass = 3
  neighbors = headers[ktmer]
  num_passed = 0
  if direction == 'forward':
    for i in range(len(neighbors)):
      if i % 3 == 2:
        if int(neighbors[i]) > dist_threshold:
          num_passed += 1
  if direction == 'backward':
    for i in range(len(neighbors)):
      if i % 3 == 1:
        if int(neighbors[i]) > dist_threshold:
          num_passed += 1
  return num_passed >= min_pass

def find_random_edge(ktmer, edges, min_dist, max_dist, ktmers):
  # Returns a list of ktmers that are within min dist, max dist
  curr_edges = edges[ktmer]
  good_edges = []
  min_edges_to_continue = 2
  for i in range(len(curr_edges)):
    if i % 2 == 1:
      if min_dist <= int(curr_edges[i]) <= max_dist or max_dist <= int(curr_edges[i]) <= min_dist:
        if curr_edges[i - 1] in ktmers and len(edges[curr_edges[i - 1]]) >= min_edges_to_continue:
          good_edges.append(curr_edges[i - 1])
  if len(good_edges) == 0:
    return ''
  else:
    return random.choice(good_edges)

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start