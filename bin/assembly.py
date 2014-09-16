# A-Bruijn Graph Construction

import sys
import string
import datetime
import random
import copy
import os
from collections import defaultdict

def main():
  reads_file = sys.argv[1]
  genome_file = sys.argv[2]
  _k = int(sys.argv[3])
  _t = int(sys.argv[4])
  gv_file = sys.argv[5]
  assembly(reads_file, genome_file, _k, _t, gv_file)

def assembly(reads, genome_file, _k, _t, gvname):
  # Given input reads in fasta format, assemble the genome
  # Procedure:
  #   1. Find all (k,T)-mers in reads
  #   2. Represent reads as list of alternating (k,T)-mers
  #   3. Construct A-Bruijn graph with reads

  print 'Assembling', reads, 'with k =', _k, 'and t =', _t
  print 'Finding kt-mers...', 
  ktmers = findKTmers(reads, _k, _t)
  print '... Done.', datetime.datetime.now()
  print 'Converting reads...'
  cReads = convertReads(reads, ktmers, _k)
  print '... Done.', datetime.datetime.now()

  # Writes distance statistics to file
  distfn = 'diststats_' + reads + '.' + str(_k) + '.' + str(_t) + '.txt'
  distfn = distfn.translate(None, '/')
  try:
    os.remove(distfn)
  except:
    pass
  with open(distfn, 'a') as f:
    for i in cReads:
      for item in i[1]:
        if item > 0:
          f.write(str(item) + '\t')

  print 'Generating graph...', datetime.datetime.now()
  headNode, allnodes = aBruijn(cReads)
  print '... Done.', datetime.datetime.now()

  print 'Joining nodes with 1 indegree/outdegree and low edge dist...'
  numJoined = 0
  numJoined = joinInOutDeg1Edges(headNode, _k)
  print '... Done.', datetime.datetime.now()


  print 'Finding many neighborhoods...'
  genome = ''
  with open(genome_file) as f:
    for i, line in enumerate(f):
      if i > 0:
        genome += line.strip()

  starting_pos = 250          # 625
  jump_length = 400           # 1250
  neighborhood_width = 500    # 1250
  neighborhood_margin = 500   # 0
  position = starting_pos
  while position < len(genome) - _k + 1:
    seedktmer = genome[position : position + _k]
    if seedktmer in allnodes:
      seednode = allnodes[seedktmer]
      neighborhood(reads, seednode, neighborhood_width, neighborhood_margin, position)
      print '... Done.', position, datetime.datetime.now()
      position += jump_length
    else:
      position += 1

  # testnode = allnodes['AAATGTTAATGGTCTGAAACGGAT']
  # neighborhood(reads, testnode, neighborhood_width, neighborhood_margin)
  print '... Done.', datetime.datetime.now()

  return


  print 'Traversing graph...'
  bulges = 0
  samektmercount = 0
  num1edge = 0
  num1edge1 = 0
  num1edge1in1 = 0
  totalkmerlen = 0
  numCombinedNodes = 0
  totalCombkmerlen = 0
  traversed = set()

  nodes = set()
  for n in headNode.outnodes:
    nodes.add(n)
  
  try:
    open(gvname, 'w').close()
  except:
    pass
  with open(gvname, 'a') as f:
    f.write('digraph G {\n')
    while len(nodes) > 0:
      # print len(nodes)
      n = nodes.pop()
      traversed.add(n)
      totalkmerlen += len(n.ktmer)
      if len(n.ktmer) > _k:
        totalCombkmerlen += len(n.ktmer)
        numCombinedNodes += 1
      # print n.ktmer, n.outDegree, n.inDegree
      if n.inDegree > 0: 
        bulges += n.inDegree - 1
      print n.ktmer + ':', n.outDegree, 'edges out and', n.inDegree, 'edges in'
      findSameKtmers = set()
      for i in range(len(n.outedges)):

        # Graphing
        f.write('\t' + n.ktmer + ' -> ' + n.outnodes[i].ktmer + ' [label = "' + str(n.outedges[i]) + '"];\n')

        if n.outnodes[i].ktmer in findSameKtmers:
          samektmercount += 1
        else:
          findSameKtmers.add(n.outnodes[i].ktmer)
        print '  dist:', n.outedges[i], n.outnodes[i].ktmer, '\treadcount:', n.outReadCount[i]
        if n.outnodes[i] not in traversed:
          nodes.add(n.outnodes[i])
      if n.outDegree == 1:
        num1edge += 1
        if n.outedges[0] == 1:
          num1edge1 += 1
        if n.outnodes[0].inDegree == 1:
          num1edge1in1 += 1
    f.write('}')

  print 'Assembled', reads, 'with k =', _k, 'and t =', _t
  print '# Nodes joined:', numJoined
  print '# Bulges:', bulges
  print '# Same edge, different dist:', samektmercount 
  print '# Nodes:', len(traversed)
  print '# Nodes w/ exactly 1 out edge:\t\t', num1edge, str(float(num1edge*100)/float(len(traversed))) + '%'
  print '# Nodes w/ exactly 1 out edge, dist 1:\t', num1edge1, str(float(num1edge1*100)/float(len(traversed))) + '%'
  print '# Edges b/w nodes of outdegree/indegree 1:\t', num1edge1in1, str(float(num1edge1in1*100)/float(len(traversed))) + '%'
  print 'Avg. kt-mer len:', str(float(totalkmerlen)/float(len(traversed)))
  
  if numCombinedNodes > 0:
    print '# Combined Nodes:', numCombinedNodes
    print 'Avg. kt-mer len for combined nodes:', str(float(totalCombkmerlen)/float(numCombinedNodes))
  else:
    print '# Combined Nodes:', numCombinedNodes
    print 'Avg. kt-mer len for combined nodes: N/A'


def neighborhood(reads, centerNode, dist, margin, position):
  # Input:
  #   node: A node in the A-bruijn graph
  #   dist: Integer
  # Output:
  #   A set of sequences corresponding to the neighborhood of dist/2 around node
  # Margin refers to searching for nodes within the neighborhood by
  #   exploring farther outside and then back in.
  # Position is the true genomic position.

  nodes = dict()        # Keys are nodes, values are position, 0 = starting pt
  traversed = set()     # Set of nodes
  collected = dict()    # Keys are nodes, values are position, 0 = starting pt
  distmin = (dist * -1)/2
  distmax = dist / 2

  current = centerNode
  currpos = 0

  while True:
    # Add all neighboring nodes if not traversed already
    for i in range(len(current.outnodes)):
      nextNode = current.outnodes[i]
      if nextNode not in traversed and nextNode not in nodes.keys():
        nodes[nextNode] = currpos + current.outedges[i]
    for i in range(len(current.innodes)):
      nextNode = current.innodes[i]
      if nextNode not in traversed and nextNode not in nodes.keys():
        nodes[nextNode] = currpos - current.inedges[i]
    
    # Remove nodes that are outside of alloted distance
    tempdict = dict()
    for key in nodes:
      if distmin - margin < nodes[key] < distmax + margin:
        tempdict[key] = nodes[key]
    nodes = tempdict

    traversed.add(current)
    if distmin < currpos < distmax:
      collected[current] = currpos

    if len(nodes) == 0:
      break
    else:
      current = nodes.keys()[0]
      currpos = nodes[current]
      del nodes[current]

    print len(nodes)

  # Find starting position and read for regions in neighborhood
  hoodReads = []
  hoodStartPos = []
  for n in collected:
    print n.ktmer, n.reads, n.pos, collected[n]
    for i in range(len(n.reads)):
      if n.reads[i] not in hoodReads:
        hoodReads.append(n.reads[i])
        hoodStartPos.append(n.pos[i] + distmin - collected[n])

  print hoodReads, hoodStartPos

  readNum = -1
  seqs = []
  with open(reads) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        readNum += 1
        if readNum in hoodReads:
          seqs.append(line.strip())
      if readNum in hoodReads and line[0] != '>':
        startpos = hoodStartPos[hoodReads.index(readNum)]
        # print readNum, line
        if startpos > 0:
          seqs.append(line[startpos:startpos + dist].strip())
        else:
          seqs.append(line[0:startpos + dist].strip())

  nhood_filename = 'e_coli_nhoods_500/nhood_nh' + str(position) + '_' + centerNode.ktmer + '.fasta'
  try:
    open(nhood_filename, 'w').close()
  except:
    pass
  with open(nhood_filename, 'a') as f:
    for seq in seqs:
      f.write(seq + '\n')

def findKTmers(reads, _k, _t):
  # Input:
  #   file: Reads in fasta format
  #   int: _k
  #   int: _t
  # Output:
  #   Set of all ktmers in reads  

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
  ans = set()

  filter_maximum = 10
  for key, val in counts.iteritems():
    if _t <= val <= filter_maximum:
      ans.add(key)
  print 'Found', len(ans), 'ktmers in', readcount, 'reads'
  return ans


def convertReads(reads, ktmers, _k):
  # Input:
  #   file: Reads in fasta format
  #   set: ktmers
  #   int: _k
  # Output:
  #   A list of tuples, each tuple is 2 lists. Each tuple describes 1 read.
  #   Each tuple contains two lists, the first list is all (k,t)-mers
  #   and the second list is the distance b/w each (k,t)-mer.
  #   The i-th element in the distance list describes the distance b/w
  #     the (i-1)th (k,t)-mer and the i-th (k,t)-mer.
  #   The 0-th element in the distance list is the distance of the first (k,t)-mer
  #     from the beginning of the string.

  isdna = False
  cReads = []
  count = 0

  with open(reads) as f:
    for i, line in enumerate(f):
      count = 0
      if isdna:
        isdna = False
        dna = line
        tempKmers = []
        tempDist = []
        for j in range(len(dna) - _k + 1):
          kmer = dna[j:j+_k]
          if kmer in ktmers:
            tempDist.append(count)
            tempKmers.append(kmer)
            count = 1
          else:
            count += 1
        cReads.append((tempKmers, tempDist))
      if line[0] == '>' or line[0] == '@':
        isdna = True

  # print cReads
  return cReads


def aBruijn(cReads):
  # Input:
  #   Reads stored as kt-mers and distances b/w kt-mers
  # Output:
  #   Builds a graph of nodes. Returns head node, an invisible node
  #   storing no information with edges to the first kt-mer in all reads.
  # This graph assembly glues nodes together simultaneously.

  head = Node('head')
  allnodes = dict()   # Keys = kmers, Values = Nodes

  # Build graph
  for j in range(len(cReads)):
    read = cReads[j]
    # print j, read
    ktmer = read[0]
    dist = read[1]
    current = head
    for i in range(len(ktmer)):
      if ktmer[i] in allnodes.keys():
        next = allnodes[ktmer[i]]
      else:
        next = Node(ktmer[i])
        allnodes[ktmer[i]] = next
      current.addOutEdge(dist[i], next)
      next.addInEdge(dist[i], current)
      next.reads.append(j)
      next.pos.append(sum(read[1][:i+1]))
      # print read[1][:i+1]
      # print next.ktmer, next.reads, next.pos
      current = next
  return head, allnodes


def joinInOutDeg1Edges(headNode, _k):
  # Joins two nodes linked by exactly 1 edge of length 1
  #
  # TODO: Currently doesn't join all nodes w/ exactly 1 edge, dist 1
  #   Find out why and fix
  # TODO: Check that node has 1 outdegree, receiving node has 1 indegree,
  #   and edge has dist 1  

  numJoined = 0
  traversed = set()
  nodes = set()
  for n in headNode.outnodes:
    nodes.add(n)
  
  while len(nodes) > 0:
    # print len(nodes)
    n = nodes.pop()
    traversed.add(n)
    # print n.ktmer, n.outDegree, n.inDegree

    # Only join nodes with indegree/outdegree of 1 and low edge distance
    if n.outDegree == 1 and n.outnodes[0].inDegree == 1 and n.outedges[0] < len(n.ktmer)/2:
      nextNode = n.outnodes[0]
      dist = n.outedges[0]
      # print dist
      # print n.ktmer[dist:], nextNode.ktmer[:len(n.ktmer)-dist]
      # Ensure that kt-mers align. Fix edges b/w nodes
      if n.ktmer[dist:] == nextNode.ktmer[:len(n.ktmer)-dist]:
        # print n.ktmer, nextNode.ktmer,
        n.ktmer += nextNode.ktmer[len(n.ktmer)-dist:]
        # print n.ktmer
        numJoined += 1
        for i in range(len(nextNode.outnodes)):
          # if nextNode.outnodes[i] not in traversed:
          #   nodes.add(nextNode.outnodes[i])
          n.addOutEdge(nextNode.outedges[i] + dist, nextNode.outnodes[i])
          nextNode.outnodes[i].addInEdge(nextNode.outedges[i] + dist, n)
          nextNode.outnodes[i].removeInEdge(nextNode.outedges[i], nextNode)
        n.removeOutEdge(n.outedges[0], nextNode)
        nextNode.removeInEdge(nextNode.inedges[0], n)
        nodes.add(n)
      else:
        print 'ERROR:', n.ktmer, nextNode.ktmer
    else:
      for i in range(len(n.outedges)):
        if n.outnodes[i] not in traversed:
          nodes.add(n.outnodes[i])

  return numJoined


class Node():
  def __init__(self, ktmer):
    self.ktmer = ktmer
    self.outReadCount = list()
    self.outedges = list()
    self.outnodes = list()
    self.inedges = list()
    self.innodes = list()
    self.outDegree = 0
    self.inDegree = 0
    self.reads = list()
    self.pos = list()

  def addOutEdge(self, edge, node):
    for i in range(len(self.outedges)):
      if node == self.outnodes[i]:
        if edge == self.outedges[i]:
          self.outReadCount[i] += 1
          return
    self.outedges.append(edge)
    self.outnodes.append(node)
    self.outReadCount.append(1)
    self.outDegree += 1

  def removeOutEdge(self, edge, node):
    if node in self.outnodes:
      for i in range(len(self.outnodes)):
        if self.outnodes[i] == node and self.outedges[i] == edge:
          del self.outnodes[i]
          del self.outedges[i]
          self.outDegree -= 1
          return

  def getAllOutEdges(self):
    return (self.outnodes, self.outedges)

  def addInEdge(self, edge, node):
    # Do not increase indegree for head
    if node.ktmer == 'head':
      return

    for i in range(len(self.inedges)):
      if node == self.innodes[i]:
        if edge == self.inedges[i]:
          return    
    self.inedges.append(edge)
    self.innodes.append(node)
    self.inDegree += 1

  def removeInEdge(self, edge, node):
    if node in self.innodes:
      for i in range(len(self.innodes)):
        if self.innodes[i] == node and self.inedges[i] == edge:
          del self.innodes[i]
          del self.inedges[i]
          self.inDegree -= 1
          return

  def getAllInEdges(self, node):
    return (self.innodes, self.inedges)

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start