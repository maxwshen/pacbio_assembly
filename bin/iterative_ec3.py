# Performs iterative error correction and production of consensus sequences in contigs.
# Beginning with an arbitrary ktmer and new contig list,
#   1. Run error correction on this ktmer
#   2. Add consensus to the contig list
#   3. Find a new ktmer some distance away from the old ktmer and repeat
#     - If there are no new ktmers, end the contig and start a new contig list at a new ktmer

# Only uses the reads in (22,4)-mers as opposed to full neighborhoods, as in iterative_ec2.py

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf
import find_read
from collections import defaultdict

def main():
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  ktmer_edges_file = '/home/mshen/research/data/22.4_ktmer_edges.out'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  min_dist = 20

  print 'Reads File:', reads_file, '\nktmer Edges File:', ktmer_edges_file, '\nktmer Headers File:', ktmer_headers_file 

  # Generate 100 sample read sets for Yu
  # with open(ktmer_headers_file) as f:
  #   lines = f.readlines()
  # for i in range(100):
  #   out_file = 'abruijn_22.4_' + str(i) + '.fasta'
  #   ktmer = lines[i].split()[0]      
  #   put_ktmer_reads_into_file(ktmer, ktmer_headers_file, reads_file, out_file)
  # return

  iterative_ec(reads_file, ktmer_headers_file, ktmer_edges_file, min_dist)

def put_ktmer_reads_into_file(ktmer, ktmer_headers_file, reads_file, out_file):
  headers = []
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      if line.split()[0] == ktmer:
        headers = line.split()[1:]
  with open(out_file, 'w+') as f:
    for h in headers:
      f.write(find_read.find_read(h, reads_file))
  return

def combine_contigs():
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'
  pass

def iterative_ec(reads_file, ktmer_headers_file, ktmer_edges_file, min_dist):
  ec_tool = '/home/mshen/research/bin/error_correction.sh'

  ktmers = []
  with open(ktmer_edges_file) as f:
    for i, line in enumerate(f):
      ktmers.append(line.split()[0])
  random.shuffle(ktmers)

  edges = build_edges_dict(ktmer_edges_file)

  restart_contig = False
  contigs = []    # List of lists
  curr_contig = []
  curr = ktmers[0]
  forward = True
  orig = curr
  while len(ktmers) > 0:
    if restart_contig:
      print len(contigs), len(ktmers)
      restart_contig = False
      contigs.append(curr_contig)
      curr_contig = []
      curr = ktmers[0]
      forward = True
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
    if forward:
      curr_contig.append(ec_seq)
    elif ec_seq not in curr_contig:   # Reinserts twice when we switch from forward to back
      curr_contig.insert(0, ec_seq)

    max_dist = calc_max_dist(curr, ec_seq, forward)

    print 'dist', max_dist, len(ec_seq)

    found = False
    edge = find_best_edge(curr, edges, min_dist, max_dist, ktmers)
    if edge != '':
      curr = edge
    else:
      if forward:
        forward = False
        curr = orig
      else:
        restart_contig = True

def calc_max_dist(curr, ec_seq, forward):
  if curr in ec_seq:
    max_dist = len(ec_seq) - ec_seq.index(curr)
    if not forward:
      max_dist -= len(ec_seq)
  else:
    max_dist = len(ec_seq) / 2
    if not forward:
      max_dist *= -1  
  return max_dist

def build_edges_dict(ktmer_edges_file):
  edges = defaultdict(list)
  with open(ktmer_edges_file) as f:
    for i, line in enumerate(f):
      edges[line.split()[0]] = line.split()[1:]
  return edges

def find_best_edge(ktmer, edges, min_dist, max_dist, ktmers):
  # Returns a list of ktmers that are within min dist, max dist
  curr_edges = edges[ktmer]
  best_edge = ''
  best_degree = 0
  for i in range(len(curr_edges)):
    print i
    if min_dist <= int(curr_edges[i]) <= max_dist or max_dist <= int(curr_edges[i]) <= min_dist:
      if len(edges[curr_edges[i - 1]]) > best_degree and curr_edges[i - 1] in ktmers:
        best_degree = len(edges[curr_edges[i - 1]])
        best_edge = curr_edges[i - 1]
  return best_edge

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start