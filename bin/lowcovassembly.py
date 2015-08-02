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
      chimerism = False
      for i, line in enumerate(f):
        words = line.split()
        if words[0] == 'Chimeric':
          self.chimeras.add(words[1])
          chimerism = True
        if words[0] == 'NonChimeric':
          chimerism = False
        if words[0] == 'Right':
          base = words[1]
          extend = words[2]
          shift = int(words[3])
          self.add_right_edge(base, extend, chimerism)

    print 'Found', len(self.chimeras), 'chimeras'
    print 'Found', len(self.nodes), 'reads'
    self.get_starting_nodes()
    self.get_ending_nodes()

    cc = self.find_connected_components()
    num_singles = 0
    for c in cc:
      sp = self.get_starting_nodes(c)
      if len(c) == 1:
        num_singles += 1
    print 'Found', num_singles, 'number of single node components'

  def add_right_edge(self, base, extend, chimerism):
    if base not in self.nodes:
      self.nodes[base] = Node(base)
    if extend not in self.nodes:
      self.nodes[extend] = Node(extend)
    self.nodes[base].add_out(extend, chimerism)
    self.nodes[extend].add_in(base, chimerism)
    return

  def get_starting_nodes(self, inp = None):
    if inp is None:
      sn = [s for s in self.nodes if len(self.nodes[s].non_inedges) == 0]
      print 'Found', len(sn), 'starting nodes in all'
    else:
      if len(inp) == 0:
        print 'Given empty list in get_starting_nodes'
        return []
      sn = [s for s in inp if len(self.nodes[s].non_inedges) == 0]
      print 'Found', len(sn), 'starting nodes in given list'
    return sn

  def get_ending_nodes(self):
    en = [s for s in self.nodes if len(self.nodes[s].non_outedges) == 0]
    print 'Found', len(en), 'ending pts found'
    return en

  def find_connected_components(self):
    used = set()    # Nums
    cc = []   # List of starting nodes

    for sn in self.get_starting_nodes():
      curr_node = self.nodes[sn]
      if curr_node in used:
        continue
      next = curr_node.non_outedges  # list of nums
      curr_cc = [curr_node.num]
      used.add(curr_node.num)
      while len(next) != 0:
        print 'Next:', len(next)
        curr_node = self.nodes[next[0]]
        next = next[1:]
        used.add(curr_node.num)
        curr_cc.append(curr_node.num)
        next += [s for s in curr_node.non_inedges if s not in used and s not in next]
        next += [s for s in curr_node.non_outedges if s not in used and s not in next]
      cc.append(curr_cc)
      print 'Used:', len(used)
    print 'Found', len(cc), 'connected components'
    return cc



class Node():
  def __init__(self, num):
    self.num = num
    self.non_outedges = []
    self.non_inedges = []
    self.chi_outedges = []
    self.chi_inedges = []

  def add_out(self, outnum, chimerism):
    if chimerism:
      self.chi_outedges.append(outnum)
    else:
      self.non_outedges.append(outnum)

  def add_in(self, innum, chimerism):
    if chimerism:
      self.chi_inedges.append(innum)
    else:
      self.non_inedges.append(innum)


if __name__ == '__main__':
  prior = '/home/yu/mshen/pacbio_assembly/data/max_assemble_files/'

  parser = argparse.ArgumentParser(description = 'Low Coverage Assembly with Overlaying Information')
  parser.add_argument('overlap_fn', \
    help = 'Reads Overlap File', \
    type = os.path.abspath, \
    nargs = '?', \
    default = prior + 'reads.20k.cov20.more_samp1.fasta.info')
  parser.add_argument('reads_fn', \
    help = 'Reads File', \
    type = os.path.abspath, \
    nargs = '?', \
    default = prior + 'reads.20k.cov20.more_samp1.fasta.info.number.fasta')
  args = parser.parse_args()
  
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main(args)
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
