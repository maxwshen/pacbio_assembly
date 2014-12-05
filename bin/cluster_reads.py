# Clusters reads in a neighborhood

import sys, string, datetime, random, copy, os, commands

import assembly
import locAL
import read_fasta
import kmer_matching

from collections import defaultdict
from subprocess import call

def main():
  reads_file = sys.argv[1]
  # directory = sys.argv[1]

  # batch(directory)

  cluster_reads_km(reads_file)
  # cluster_reads(reads_file)
  return

def batch(directory):
  files = os.listdir(directory)
  for fil in files:
    print directory + fil, datetime.datetime.now()
    cluster_reads(directory + fil)

def cluster_reads_km(reads_file):
  h, r = read_fasta.read_fasta(reads_file)
  clusters = []

  while len(r) > 0:
    r1 = r[0]
    _k = 15
    cutoff = 5
    headers = kmer_matching(r1, h[1:], r[1:], _k, cutoff, len(clusters))
    clusters.append(headers.append(h[0]))
    for item in headers:
      del r[h.index(item)]
      del h[h.index(item)]
    print len(r)


def kmer_matching(ec_seq, hr, rr, _k, cutoff, num):
  kmers = set()
  for i in range(len(ec_seq) - _k + 1):
    kmers.add(ec_seq[i:i + _k])

  reads = dict()    # Key = header, value = num shared kmers
  for i in range(len(rr)):
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  headers = []
  for key in sorted(reads, key = reads.get, reverse = True):
    if reads[key] < cutoff:
      break
    headers.append(key)

  out_file = 'temp.fasta'
  to_write = ''
  for h in headers:
    to_write += h + '\n' + rr[hr.index(h)]
  with open(out_file, 'w') as f:
    f.write(to_write)

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  blasr_options = '-bestn 1'
  blasr_out = commands.getstatusoutput(blasr_exe + ' ' + out_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]

  to_print = []
  for line in blasr_out.splitlines():
    head = '>' + '/'.join(line.split()[0].split('/')[:-1])
    start_pos = line.split()[6]
    end_pos = line.split()[7]
    length = int(end_pos) - int(start_pos)
    info = (head, str(reads[head]), start_pos, end_pos, str(length))
    to_print.append(info)

  for t in sorted(to_print, key = lambda tup: int(tup[1]), reverse = True):
    print t[0]
    print str(num) + '\t' + '\t'.join(t[1:])

  return headers


def cluster_reads(reads_file):
  h, r = read_fasta.read_fasta(reads_file)

  len_cutoff = 400
  r = filter_by_len(r, len_cutoff)

  clusters = []

  while len(r) > 0:
    r1 = r[0]
    curr_cluster = [h[0], r1]
    for j in range(1, len(r)):
      r2 = r[j]
      (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external(r1, r2, 1, -2, -2, -1)
      # print j, alignLen
      if alignLen > len_cutoff:
        curr_cluster.append(h[j])
        curr_cluster.append(r[j])
    clusters.append(curr_cluster)
    for item in curr_cluster:
      if item[0] == '>':
        del r[h.index(item)]
        del h[h.index(item)]
    # print curr_cluster, '\n', len(curr_cluster)
    # print len(r), len(h)    

  for i in range(len(clusters)):
    out_file = 'clusters_22.4_100k/nhood_' + reads_file.split('_')[6] + '_cluster_' + str(i) + '.fasta'
    with open(out_file, 'w') as f:
      for line in clusters[i]:
        f.write(line + '\n')

  return


def filter_by_len(reads, cutoff):
  new_r = []
  for r in reads:
    if len(r) > cutoff:
      new_r.append(r)
  return new_r

def longest_common_substring(s1, s2):
   m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
   longest, x_longest = 0, 0
   for x in xrange(1, 1 + len(s1)):
       for y in xrange(1, 1 + len(s2)):
           if s1[x - 1] == s2[y - 1]:
               m[x][y] = m[x - 1][y - 1] + 1
               if m[x][y] > longest:
                   longest = m[x][y]
                   x_longest = x
           else:
               m[x][y] = 0
   return s1[x_longest - longest: x_longest]

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start