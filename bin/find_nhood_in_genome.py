# Uses blasr to find where nhoods are in the genome

import sys
import string
import datetime
import random
import copy
import os
import commands
import numpy as np

def main():
  # read_file = sys.argv[1]
  directory = sys.argv[1]
  genome_file = sys.argv[2]

  batch(directory, genome_file)
  # find_nhood_in_genome(read_file, genome_file)
  return

def batch(directory, genome_file):
  clusters_per_read_file = 'clusters_per_read.out'
  cluster_sizes_file = 'cluster_sizes.out'

  files = os.listdir(directory)
  for fil in files:
    clusters = find_nhood_in_genome(fil, genome_file)
    with open(cluster_sizes_file, 'a') as f:
      for c in clusters:
        f.write(str(c))
    with open(clusters_per_read_file, 'a') as f:
      f.write(str(len(c)))


def find_nhood_in_genome(read_file, genome_file):
  region_width = 1000
  blasr = '/home/jeyuan/blasr/alignment/bin/blasr'
  sam = commands.getstatusoutput(blasr + ' ' + read_file + ' ' + genome_file + ' -bestn 1')[1]
  pos = []
  for line in sam.splitlines():
    # print line
    pos.append(int(line.split()[6]))

  clusters = []   # List of lists
  curr = []
  while len(pos) > 0:
    curr = [pos[0]]
    for j in range(1, len(pos)):
      if np.absolute(pos[0] - pos[j]) < region_width:
        curr.append(pos[j])
    for c in curr:
      del pos[pos.index(c)]
    clusters.append(curr)

  return [len(c) for c in clusters]


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start