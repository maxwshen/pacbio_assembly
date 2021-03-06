# 

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict

def main(args):
  overlap_fn = args.overlap_fn
  reads_fn = args.reads_fn

  og = OverlapGraph(overlap_fn)
  og.find_nonchimeric_contigs()
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
          right_shift = int(words[3])
          self.add_right_edge(base, extend, chimerism, right_shift)

    print 'Found', len(self.chimeras), 'chimeras'
    print 'Found', len(self.nodes), 'reads'
    self.get_starting_nodes()
    self.get_ending_nodes()

  def find_nonchimeric_contigs(self):
    contigs = []
    cc = self.find_connected_components()
    num_singles = 0
    nums = dict()
    for c in cc:
      if len(c) == 1:
        num_singles += 1
        contigs.append(c)
        continue
      sp = self.get_starting_nodes(c)
      self.longest_path(sp)

    print '...Found', num_singles, 'single node components out of', len(cc)
    for key in nums:
      print '...Found', nums[key], 'components with', key, 'starting nodes'
    return

  def longest_path(self, starting_nodes):
    traversed = set()   # Stores nums
    best = dict()   # Key = node num, val = longest distance
    queue = starting_nodes
    while len(queue) != 0:
      curr = self.nodes[queue[0]]
      queue = queue[1:]
      traversed.add(curr.num)

      # Update scores and add ready nodes to queue
      if curr.num not in best:
        base = 0    # Technically shouldn't be 0, should be length of read
      else:
        base = best[curr.num]
      for oe in curr.non_outedges:
        if oe in traversed or oe in queue:
          continue
        new_score = base + curr.outweights[oe]
        if oe not in best:
          best[oe] = new_score
        elif new_score > best[oe]:
          best[oe] = new_score
        if self.nodes[oe].is_ready(traversed):
          if oe in queue:
            queue.append(oe)

    print 'Num traversed:', len(traversed)
    ends = self.get_ending_nodes(list(traversed))
    print 'Num nodes in best:', len(best), 
    print 'Num ends:', len(ends) 
    if len(ends) == 0:
      print 'No ends found - circle?'
    else:
      print max([best[s] for s in ends]), [best[s] for s in ends]
    return

  def add_right_edge(self, base, extend, chimerism, right_shift):
    if base not in self.nodes:
      self.nodes[base] = Node(base)
    if extend not in self.nodes:
      self.nodes[extend] = Node(extend)
    self.nodes[base].add_out(extend, chimerism, right_shift)
    self.nodes[extend].add_in(base, chimerism)
    return

  def get_starting_nodes(self, inp = None):
    # Inp can be a list of numbers
    if inp is None:
      sn = [s for s in self.nodes if len(self.nodes[s].non_inedges) == 0]
      print 'Found', len(sn), 'starting nodes in all'
    else:
      if len(inp) == 0:
        print 'Given empty list in get_starting_nodes'
        return []
      sn = [s for s in inp if len(self.nodes[s].non_inedges) == 0]
      # print 'Found', len(sn), 'starting nodes in given list'
    return sn

  def get_ending_nodes(self, inp = None):
    # Inp can be a list of numbers
    if inp is None:
      en = [s for s in self.nodes if len(self.nodes[s].non_outedges) == 0]
      print 'Found', len(en), 'ending pts found in all'
    else:
      if len(inp) == 0:
        print 'Given empty list in get_starting_nodes'
        return []
      en = [s for s in inp if len(self.nodes[s].non_outedges) == 0]
    return en

  def find_connected_components(self):
    used = set()    # Nums
    cc = []   # List of lists, each list contains node nums

    for sn in self.get_starting_nodes():
      curr_node = self.nodes[sn]
      if curr_node in used:
        continue
      next = curr_node.non_outedges  # list of nums
      curr_cc = [curr_node.num]
      used.add(curr_node.num)
      while len(next) != 0:
        # print 'Next:', len(next)
        curr_node = self.nodes[next[0]]
        next = next[1:]
        used.add(curr_node.num)
        curr_cc.append(curr_node.num)
        next += [s for s in curr_node.non_inedges if s not in used and s not in next]
        next += [s for s in curr_node.non_outedges if s not in used and s not in next]
      cc.append(curr_cc)
      # print 'Used:', len(used)
    print 'Found', len(cc), 'connected components'
    return cc



class Node():
  def __init__(self, num):
    self.num = num
    self.non_outedges = []
    self.outweights = dict()   # Key = neighbor, val = right shift int
    self.non_inedges = []
    self.chi_outedges = []
    self.chi_inedges = []

  def is_ready(self, traversed):
    return len(traversed.intersection(self.non_inedges)) == len(traversed)

  def add_out(self, outnum, chimerism, right_shift):
    if chimerism:
      self.chi_outedges.append(outnum)
    else:
      self.non_outedges.append(outnum)
    if outnum in self.outweights:
      print 'ERROR:', outnum, 'already in self.outweights of node', self.num
    else:
      self.outweights[outnum] = right_shift

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
