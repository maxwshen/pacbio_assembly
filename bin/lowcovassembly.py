# 

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict

def main(args):
  overlap_fn = args.overlap_fn
  reads_fn = args.reads_fn

  og = OverlapGraph(overlap_fn)
  return


class OverlapGraph():
  def __init__(self, overlap_fn):
    self.chimeras = set()
    self.nodes = dict()   # Key = num, Val = node

    with open(overlap_fn) as f:
      for i, line in enumerate(f):
        words = line.split()
        if words[0] == 'Chimeric':
          self.chimeras.add(words[1])
        if words[0] == 'Right':
          base = words[1]
          extend = words[2]
          shift = int(words[3])

          self.add_right_edge(base, extend)

    print 'Found', len(self.chimeras), 'chimeras'
    print 'Found', len(self.nodes), 'reads'
    print len([s for s in self.nodes if len(self.nodes[s].inedges) == 0]), 'starting pts found'

  def add_right_edge(self, base, extend):
    if base not in self.nodes:
      self.nodes[base] = Node(base)
    if extend not in self.nodes:
      self.nodes[extend] = Node(extend)
    self.nodes[base].add_out(extend)
    self.nodes[extend].add_in(base)
    return

class Node():
  def __init__(self, num):
    self.num = num
    self.outedges = []
    self.inedges = []

  def add_out(self, outnum):
    self.outedges.append(outnum)

  def add_in(self, innum):
    self.inedges.append(innum)



if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = 'Low Coverage Assembly with Overlaying Information')
  parser.add_argument('overlap_fn', \
    help = 'Reads Overlap File', \
    type = os.path.abspath, \
    nargs = '?', \
    default = '/home/yu/mshen/pacbio_assembly/data/max_assemble_files/reads.20k.cov20.more_samp1.fasta.info')
  parser.add_argument('reads_fn', \
    help = 'Reads File', \
    type = os.path.abspath, \
    nargs = '?', \
    default = '/home/yu/mshen/pacbio_assembly/data/max_assemble_files/reads.20k.cov20.more_samp1.fasta.info.number.fasta')
  args = parser.parse_args()
  
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main(args)
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start