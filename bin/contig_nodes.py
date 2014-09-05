# contig.py
# Modified from CSE 181
# 
# Relies on assembly.py and locAL.py
#
# Finds <k> <t> mers in an input fasta file and forms contigs 
# using a debruijn graph assembly approach.
#
# Input file and genome expected to be in fasta format.
# Genome expected to just have 1 sequence
#
# Uses a node class as opposed to an adjacency list

import sys
import string
import datetime
import random
import copy
import assembly
import locAL
import os

from collections import defaultdict

def main():
  global kparam
  global tparam
  global lparam
  global trimparam
  global silenced
  silenced = False

  try:
    reads_file = sys.argv[1]
    genome_file = sys.argv[2]
    kparam = int(sys.argv[3])
    tparam = int(sys.argv[4])
    gvname = sys.argv[5]
    lparam = int(sys.argv[6])
    if lparam >= kparam:
      print 'ERROR: Bad l: Greater than or equal to k'
      sys.exit()
    trimparam = int(sys.argv[7])
    if trimparam < 1:
      trimparam = 0
  except:
    print 'Usage: python contig_nodes <input file> <genome file> <k> <t> <graphviz name> <l param> <trimmingNum>'
    sys.exit(0)


  use_high_deg_nodes = True
  if use_high_deg_nodes:
    kmers = []
    with open(reads_file) as f:
      for i, line in enumerate(f):
        kmers.append((line.split()[0], int(line.split()[3]), float(line.split()[10])))
  else:
    kmers = findKTmers(reads_file, kparam, tparam)

  if lparam != 0:
    kmers = splitKTmers(kmers, kparam - lparam)
  
  nodes = deBruijnkmers(kmers)
  nodes = trimdebruijn(nodes, trimparam)
  contigs = contig(nodes, genome_file, gvname, kparam, lparam)
  supercontigs = findPathsNew(nodes, kparam)
  find_location(supercontigs, nodes, kparam)
  checkAccuracy(supercontigs, genome_file)

  # bestpath = findPaths(nodes)
  # checkAccuracy([bestpath], genome_file)

def find_location(super_contigs, nodes, _k):
  positional_multiplier = 0.916
  klen = _k - 1
  for contig in super_contigs:
    start = contig[:klen]
    end = contig[-klen:]
    startpos = nodes[start].pos * positional_multiplier
    endpos = nodes[end].pos * positional_multiplier + _k
    print '[' + str(startpos) + '] [' + str(endpos) + '] ' + contig

  for contig in super_contigs:
    sum_deg = 0
    for i in range(len(contig) - klen + 1):
      kmer = contig[i : i + klen]
      sum_deg += nodes[kmer].tval
    print contig, float(sum_deg) / (len(contig) - klen)



def test_canonical_paths(reads, genome_file, _k, _t, gvname, _l, _trim, ktmers, silenceFlag = True):
  # Used to call contig from another script.
  # ktmers is a set of strings, not necessarily same k, t as _k, _t
  global silenced
  silenced = silenceFlag

  if not os.path.isfile(reads) or not os.path.isfile(genome_file):
    print 'ERROR: Bad reads or genome file'
    return 0

  global kparam
  global tparam
  global lparam
  global trimparam
  kparam = _k
  tparam = _t
  lparam = _l
  trimparam = _trim

  kmers = findKTmers(reads, kparam, tparam)
  if lparam != 0:
    kmers = splitKTmers(kmers, kparam - lparam)

  nodes = deBruijnkmers(kmers)

  nodes = trimdebruijn(nodes, trimparam)
  contigs = contig(nodes, genome_file, gvname)
  # bestpath = findPaths(nodes)
  # score, alignLen = checkAccuracy([bestpath], genome_file)  

  # return score, alignLen

  trimmed_ktmers = []
  missing = []
  for ktmer in ktmers:
    if findKmer(ktmer, nodes):
      trimmed_ktmers.append(ktmer[:_k - 1])
    else:
      missing.append(ktmer)

  num_true = 0
  num_false = 0

  for i in range(len(trimmed_ktmers)):
    for j in range(i, len(trimmed_ktmers)):
      print trimmed_ktmers[i], trimmed_ktmers[j],
      if isConnected(nodes, trimmed_ktmers[i], trimmed_ktmers[j]):
        print 'True'
        num_true += 1
      else:
        print 'False'
        num_false += 1

  if num_true + num_false > 0:
    print float(num_true) / float(num_true + num_false)
  else:
    print 'None found'

  # print 'Missing:', missing

def isConnected(nodes, node1kmer, node2kmer):
  # Returns true/false if node1 and node2 are connected
  curr = nodes[node1kmer]
  traversed = set()
  pending = []
  for n in curr.inedges:
    pending.append(n)
  for n in curr.outedges:
    pending.append(n)
  traversed.add(curr)

  while len(pending) > 0:
    curr = pending[0]
    pending = pending[1:]
    traversed.add(curr)
    for n in curr.inedges:
      if n.kmer == node2kmer:
        return True
      if n not in traversed:
        pending.append(n)
    for n in curr.outedges:
      if n.kmer == node2kmer:
        return True
      if n not in traversed:
        pending.append(n)

  return False

def findKmer(text, nodes):
  # Given a text dna word, finds it in the nodes if it exists
  # Returns true or false

  klen = len(nodes.keys()[0])
  text_start = text[:klen]

  try:
    curr = nodes[text_start]
  except:
    return False

  for char in text[klen:]:
    found = False
    for neighbor in curr.outedges:
      if neighbor.kmer[-1] == char:
        curr = neighbor
        found = True
    if not found:
      return False

  return True

def contig_nodes(reads, genome_file, _k, _t, gvname, _l, _trim, silenceFlag = True):
  # Used to call contig from another script
  global silenced
  silenced = silenceFlag

  if not os.path.isfile(reads) or not os.path.isfile(genome_file):
    print 'ERROR: Bad reads or genome file'
    return 0

  global kparam
  global tparam
  global lparam
  global trimparam
  kparam = _k
  tparam = _t
  lparam = _l
  trimparam = _trim

  kmers = findKTmers(reads, kparam, tparam)
  if lparam != 0:
    kmers = splitKTmers(kmers, kparam - lparam)

  nodes = deBruijnkmers(kmers)

  nodes = trimdebruijn(nodes, trimparam)
  contigs = contig(nodes, genome_file, gvname, kparam, lparam)
  # bestpath = findPaths(nodes)
  # score, alignLen = checkAccuracy([bestpath], genome_file)  

  # return score, alignLen

def diameter(nodes):
  # Finds the approximate diameter of the graph
  starting_nodes = []
  for n in nodes.values():
    if n.outdegree == 1 and n.indegree == 0:
      starting_nodes.append(n)

  diameter = 0
  for n in starting_nodes:
    temp_diameter = 1
    curr = n
    while len(curr.outedges) != 0:
      curr = curr.outedges[0]
      temp_diameter += 1
    if temp_diameter > diameter:
      diameter = temp_diameter

  return diameter

def findPathsNew(nodes, _k):
  # Given a dict of nodes (key = kmer, value = Node object), finds the longest path
  # for each disconnected component of the de Bruijn graph. Basically it's like a super-contig,
  # a way to trim erroneous tips.
  global silenced
  starting_nodes = []
  for n in nodes.values():
    if n.outdegree == 1 and n.indegree == 0:
      starting_nodes.append(n)

  bestpaths = []
  for n in starting_nodes:
    paths = [[n]]   # A list of lists of nodes. Most recent node in a path is at end
    traversed = dict()    # Key = Node, Value = Lists containing the node
    completed_paths = 0

    while True:
      completed_paths = 0
      for path in paths:
        if path[-1] == 'End':
          completed_paths += 1
      if completed_paths == len(paths):
        # print [p.kmer for p in path if not isinstance(p, basestring)]
        bestpath = paths[0]
        for path in paths:
          if len(path) > len(bestpath):
            bestpath = path
        bestpaths.append(bestpath)
        break

      new_paths = []
      for i in range(len(paths)):
        path = paths[i]
        if path[-1] != 'End':
          curr = path[-1]  
          if len(curr.outedges) == 1:
            neighbor = curr.outedges[0]
            if neighbor not in path:
              path.append(neighbor)
              if neighbor not in traversed:
                traversed[neighbor] = i
              else:
                paths[i], paths[traversed[neighbor]] = mergeBulge(paths[i], paths[traversed[neighbor]])
            else:
              path.append('End')    # Encountered a loop
          elif len(curr.outedges) == 0:
            path.append('End')
          elif len(curr.outedges) > 1:
            for neighbor in curr.outedges:
              temp_path = copy.copy(path)
              temp_path.append(neighbor)
              new_paths.append(temp_path)
            path.append('End')

      for new_path in new_paths:
        paths.append(new_path)

  print 'Done'

  path_strings = []
  for path in bestpaths:
    path_string = path[0].kmer
    for i in range(1, len(path) - 1):
      path_string += path[i].kmer[-1]
    path_strings.append(path_string)
    for text in path_strings:
      if text[-_k:] == path_string[-_k:] and text != path_string:
        if len(path_string) > len(text):
          path_strings.remove(text)
        else:
          path_strings.remove(path_string)

  for p in path_strings:
    print p

  return path_strings

def mergeBulge(list1, list2):
  # Input: Two lists of nodes that have two nodes in common
  # Puts an 'End' on one of the lists and returns them.
  # Calculates the average t of the nodes in the bulge

  if list1[-1] == 'End':
    list1 = list1[:-1]
  if list2[-1] == 'End':
    list2 = list2[:-1]

  end_node = list1[-1]
  list1_end_index = len(list1) - 1
  list2_end_index = list2.index(end_node)

  print 'Bulge:', [n.kmer for n in list1 if not isinstance(n, basestring)]
  print 'Bulge:', [n.kmer for n in list2 if not isinstance(n, basestring)]

  for i in range(len(list1)):
    for j in range(len(list2)):
      if list1[i].kmer == list2[j].kmer and j != list2_end_index:
        list1_start_index = i
        list2_start_index = j
        start_node = list1[i]

  print list1_start_index, list1_end_index, list2_start_index, list2_end_index
  print list1[list1_start_index].kmer, list1[list1_end_index].kmer
  print list2[list2_start_index].kmer, list2[list2_end_index].kmer

  sum_t1 = 0
  sum_t2 = 0
  num_t1 = 0
  num_t2 = 0
  avg_t1 = 0
  avg_t2 = 0

  for i in range(list1_start_index, list1_end_index):
    print list1[i].kmer, list1[i].tval
    sum_t1 += int(list1[i].tval)
    num_t1 += 1
  for i in range(list2_start_index, list2_end_index):
    print list2[i].kmer, list2[i].tval
    sum_t2 += int(list2[i].tval)
    num_t2 += 1
  avg_t1 = float(sum_t1) / float(num_t1)
  avg_t2 = float(sum_t2) / float(num_t2)

  if avg_t1 >= avg_t2:
    newlist1 = list1
    newlist2 = copy.copy(list2)
    newlist2.append('End')
  else:
    newlist1 = copy.copy(list1)
    newlist1.append('End')
    newlist2 = list2

  print avg_t1, avg_t2

  return newlist1, newlist2

def findPaths(nodes):
  # Given a dict of nodes (key = kmer, value = Node object), find longest path
  global silenced
  starting_nodes = []
  for n in nodes.values():
    if n.outdegree == 1 and n.indegree == 0:
      starting_nodes.append(n)

  bestpath = ''
  for n in starting_nodes:
    path = n.kmer
    curr = n
    while len(curr.outedges) != 0:
      if len(curr.outedges) > 1:
        text, new_node = exploreBulgeBranch(curr)
        path += text
        curr = new_node
      elif len(curr.outedges) == 1:
        curr = curr.outedges[0]
        path += curr.kmer[-1]

    if len(path) > len(bestpath):
      bestpath = path

  if not silenced:
    print 'best:', bestpath
  return bestpath

def exploreBulgeBranch(node):
  # Given a node with more than 1 outedge, explore to see if it's a branch or a bulge
  # If it's a branch, return the final node of the longer path and the text
  # If it's a bulge, return the final node joining the bulge and the text of the path
  # with higher t
  global silenced

  children = copy.copy(node.outedges)
  traversed = set()             # Set of nodes

  totalt = []
  texts = []
  paths = []    # List of lists, index = when node was traversed
  numNodes = []
  for child in children:
    texts.append(child.kmer[-1])
    totalt.append(child.tval)
    traversed.add(child)
    numNodes.append(1)
    paths.append([child.kmer])

  if not silenced:
    print node.kmer

  bulge = False
  branch = False
  while True:
    badbranches = []
    for i in range(len(children)):
      curr = children[i]
      if len(curr.outedges) == 1:
        if curr.outedges[0] not in traversed:
          curr = curr.outedges[0]
          texts[i] += curr.kmer[-1]
          totalt[i] += curr.tval
          traversed.add(curr)
          paths[i].append(curr.kmer)
          children[i] = curr
          numNodes[i] += 1
        else:
          # for t in traversed:
          #   print 't:', t.kmer
          if not silenced:
            print 'already traversed:', curr.outedges[0].kmer
          bulge = True
          bulgept = curr.outedges[0]
          for j in range(len(numNodes)):
            if j != i:
              numNodes[j] = paths[j].index(curr.outedges[0].kmer)
      elif len(curr.outedges) > 1:
        if not silenced:
          print 'WARNING:', curr.kmer, 'is branch/bulge within a branch/bulge'
        pass
        # Recursive case of exploreBulgeBranch
        # Right now, any path that represents a bulge/branch within a bulge is not
        # further explored. children[] is not updated.
        # THIS WILL ENTER AN INFINITE LOOP if all paths inside a bulge are also bulges
      elif len(curr.outedges) == 0:
        badbranches.append(i)


      if bulge:   # Don't break on branch so full iteration can be completed, all nodes extended same amt
        break

    for i in badbranches:
      del children[i]
      del texts[i]
      del totalt[i]
      del numNodes[i]
      if len(children) == 1:
        branch = True

    if branch or bulge:
      break

  if bulge:
    avgt = []
    for i in range(len(totalt)):
      avgt.append(totalt[i] / numNodes[i])
    bestindex = avgt.index(max(avgt))
    if not silenced:
      print avgt, totalt, numNodes
      print 'Best Bulge path:', texts[bestindex][:numNodes[bestindex] + 1], bulgept.kmer
    return texts[avgt.index(max(avgt))], bulgept

  if branch:
    if not silenced:
      print 'Best branch path:', texts[0], children[0].kmer
    return texts[0], children[0]

def trimdebruijn(nodes, num):
  # Input:
  #   Dictionary of nodes
  # Output:
  #   Trimmed dictionary of nodes
  # Procedure:
  #   Remove all nodes with zero indegree or outdegree, iterate as desired
  global silenced

  for i in range(num):
    kmers_to_delete = set()
    for kmer in nodes.keys():
      n = nodes[kmer]
      if n.indegree == 0:
        kmers_to_delete.add(kmer)
    for kmer in nodes.keys():
      n = nodes[kmer]
      if n.outdegree == 0:
        kmers_to_delete.add(kmer)

    for kmer in kmers_to_delete:
      n = nodes[kmer]
      if n.indegree == 0:
        for m in n.outedges:
          m.removeInEdge(n)
      elif n.outdegree == 0:
        for m in n.inedges:
          m.removeOutEdge(n)
      del nodes[kmer]

    trimcount = len(kmers_to_delete)
    if not silenced:
      print 'Trimmed', trimcount, 'leaves. New graph size:', len(nodes), 'nodes'

  return nodes

def findKTmers(reads, _k, _t):
  # Input:
  #   file: Reads in fasta format
  #   int: _k
  #   int: _t
  # Output:
  #   Set of all ktmers in reads in tuples, (dna, t)  

  global silenced
  isdna = False
  counts = dict()
  readcount = 0
  with open(reads) as f:
    for i, line in enumerate(f):
      # print i
      if isdna:
        isdna = False
        dna = line.strip()
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j+_k]
          # print kmer
          if kmer in counts:
            counts[kmer] = counts[kmer] + 1
          else:
            counts[kmer] = 1
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        isdna = True
  ans = []

  for key, val in counts.iteritems():
    if val >= _t:
      ans.append((key, val))
  if not silenced:
    print 'Found ' + str(len(ans)) + ' (' + str(_k) + ',' + str(_t) + ')-mers in ' + str(readcount) + ' reads'
  return ans

def contig(nodes, genome_file, gvname, kparam, lparam):
  global silenced

  # Visualize graph. Black edges = true, red edges = erroneous
  genome_kmers = set()
  with open(genome_file) as f:
    lines = f.readlines()
    genome = lines[1]
  for i in range(len(genome) - kparam + lparam + 1):
    genome_kmers.add(genome[i : i + kparam - lparam])

  try:
    open(gvname, 'w').close()
  except:
    pass
  with open(gvname, 'a') as f:
    f.write('digraph G {\n')
    for n in nodes.values():
      for neighbor in n.outedges:
        f.write('"' + n.kmer + '\\n' + str(n.pos) + '" -> ')
        f.write('"' + neighbor.kmer + '\\n' + str(neighbor.pos) + '"')
        if n.kmer + neighbor.kmer[-1] in genome_kmers:
          f.write(';\n')
        else:
          f.write(' [color = "red"];\n')        
    f.write('}')

  # Get contigs
  separated_parts = []
  contigs = []
  contig = ''
  for n in nodes.values():
    if n.outdegree == 1 and n.indegree == 0:
      separated_parts.append(n)

  starting_nodes = []
  for n in separated_parts:
    starting_nodes.append(n)

  # starting_nodes: nodes with outdegree 1 with no indegree or indegree a branchpoint.
  # traversed: all traversed nodes. We ensure that we do not traverse them again. 
  # contig: The string contig, built by adding the last character of each subsequent node in
  #   a non-branching path. Contig is appended to the list of contigs when we hit a branchpoint
  # old_current: Used to detect if all neighbors are traversed. If so, then we append the current
  #   contig to the list of contigs and quit that path
  traversed = set()
  while len(starting_nodes) != 0:
    # print len(starting_nodes)
    current = starting_nodes[0]
    starting_nodes = starting_nodes[1:]
    contig = current.kmer
    while len(current.outedges) != 0:
      traversed.add(current)
      old_current = copy.copy(current)
      for neighbor in current.outedges:
        if neighbor not in traversed:
          current = neighbor
      if current.kmer != old_current.kmer:
        if len(current.outedges) == 1 and len(current.inedges) == 1:
          contig += current.kmer[-1]
        elif len(current.outedges) > 1 or len(current.inedges) > 1:
          starting_nodes.append(current)
          contigs.append(contig)
          contig = ''
          break
        elif len(current.outedges) == 0:
          contig += current.kmer[-1]
          contigs.append(contig)
          break
      else:
        if len(contig) > 0:
          contigs.append(contig)
        break

  if not silenced:
    for c in contigs:
      print c

  return contigs

def checkAccuracy(seqs, genome_file):
  # Check the accuracy of the best path
  global silenced
  if not silenced:
    print 'Aligning to genome...'
  with open(genome_file) as f:
    lines = f.readlines()
    genome = lines[1]

  numPerfect = 0
  matchScores = []
  perfectStarts = []
  perfectLens = []
  for i in range(len(seqs)):
    # print seqs[i]
    (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(seqs[i], genome, 1, -1, -1, -0.5)
    if alignLen != 0:
      score = float(matches) / float(alignLen)
    else:
      score = 0
      print 'ERROR: Divide by Zero! Alignment Length = 0'
    if not silenced:
      print 'Range:', bestxy[1] - alignLen, '-', bestxy[1]
      print 'Accuracy:', score
      print 'Error Rate:', 1 - score
    matchScores.append((score, alignLen))

  return score, alignLen

def deBruijnkmers(kmers_in):
  # Input:
  #   A list of kmers (kmer, _t, pos)
  # Output:
  #   Dictionary. Keys = (k-1)mers, Values = Node object

  existing_kmers = set()
  nodes = dict() # Key = kmer, Value = Node object
  for (kmer, _t, pos) in kmers_in:
    prefix = kmer[:len(kmer)-1]
    suffix = kmer[1:]
    if prefix not in nodes.keys():
      curr_node = Node(prefix, _t, pos)
      nodes[prefix] = curr_node
    else:
      nodes[prefix].addT(_t)
    if suffix not in nodes.keys():
      next_node = Node(suffix, _t, pos)
      nodes[suffix] = next_node
    else:
      nodes[suffix].addT(_t)

    nodes[prefix].addOutEdge(nodes[suffix])
    nodes[suffix].addInEdge(nodes[prefix])

  return nodes

def overlaps(q1, q2):
  # Returns true if suffix of q1 overlaps with the prefix of q2
  # Used in deBruijnkmers
  q1 = q1[1:]
  q2 = q2[:len(q2)-1]
  if len(q1) != len(q2):
    return False
  for i in range(len(q1)):
    if q1[i] != q2[i]:
      return False
  return True

def splitKTmers(kmers, l):
  # Input:
  #   Set of kmers
  # Output:
  #   Set of l-mers created from the kmers 
  print 'Splitting ktmers with l =', l
  newkmers = []
  for i in range(len(kmers)):
    for j in range(len(kmers[i][0]) - l + 1):
      newkmers.append((kmers[i][0][j:j+l], kmers[i][1]))

  # print newkmers
  return newkmers


class Node():
  def __init__(self, kmer, tval, pos):
    self.kmer = kmer
    self.tval = tval
    self.pos = pos
    self.outedges = list()
    self.inedges = list()
    self.degree = 0     # Out edges are positive, in edges are negative
    self.outdegree = 0
    self.indegree = 0

  def addOutEdge(self, edge):
    if edge in self.outedges:
      return
    self.outedges.append(edge)
    self.degree += 1
    self.outdegree += 1

  def removeOutEdge(self, edge):
    if edge in self.outedges:
      del self.outedges[self.outedges.index(edge)]
      self.degree -= 1
      self.outdegree -= 1

  def addT(self, _t):
    self.tval += _t

  def getAllOutEdges(self):
    return self.outedges

  def addInEdge(self, edge):
    if edge in self.inedges:
      return
    self.inedges.append(edge)
    self.degree -= 1
    self.indegree += 1

  def removeInEdge(self, edge):
    if edge in self.inedges:
      del self.inedges[self.inedges.index(edge)]
      self.degree += 1
      self.indegree -= 1

  def getAllInEdges(self, node):
    return self.inedges


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start